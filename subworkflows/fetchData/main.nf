// Subworkflow to fetch sample and reference data from --fasta/--reads/--ref_fasta/--ref_reads

// Set path variables
output_directory = file(params.output_directory)
assembly_directory = file(params.assembly_directory)
log_directory = file(params.log_directory)

ref_id_file = file(params.ref_id_file)

// Set ref_mode
ref_mode = params.ref_mode

// Set file headers
assembly_header = "Isolate_ID\tRead_Type\tRead_Location\tAssembly_Path\n"

// Set paths to accessory files/scripts
assembly_log = file("${log_directory}/Assembly_Data.tsv")
user_snpdiffs_list = file("${log_directory}/Imported_SNPDiffs.txt")
findReads = file("${projectDir}/bin/fetchReads.py")
userSNPDiffs = file("${projectDir}/bin/userSNPDiffs.py")

// Set SKESA cores to 5 or fewer
skesa_cpus = (params.cores as Integer) >= 5 ? 5 : (params.cores as Integer)

workflow {
    main:
    query_data, reference_data, snpdiff_data = fetchData()

    publish:
    query_data >> 'query_data.tsv'
    reference_data >> 'reference_data.tsv'
    snpdiff_data >> 'snpdiff_data.tsv'
}

// Top-level workflow //
workflow fetchData{

    emit:
    query_data
    reference_data
    snpdiff_data

    main:
    // Get any excluded IDs
    ("${params.exclude}" != "" ? processExclude() : Channel.empty()).set{exclude_ids} 
    
    // Process snpdiffs alignments
    // If assembly file cannot be found, it will be 'null'
    ("${params.snpdiffs}" != "" ? processSNPDiffs() : Channel.empty()).set{user_snpdiffs}

    excluded_snpdiffs = user_snpdiffs.map{it -> tuple(it[1],it[0])}
    .concat(user_snpdiffs.map{it -> tuple(it[10],it[0])})
    .join(exclude_ids,by:0)
    .unique{it -> it[1]}
    .map{it -> tuple(it[1],"Exclude")}

    // Generate return channel: 3-item tuple (Query_ID, Reference_ID, SNPDiff_Path)
    snpdiff_data = user_snpdiffs
    .map{it -> tuple(it[0],it[1],it[10])}
    .join(excluded_snpdiffs,by:0,remainder:true)
    .filter{it -> it[0].toString() != "null"}
    .filter{it -> it[3].toString() != "Exclude"}
    .unique{it -> it[0]}
    .map{it -> tuple(it[1],it[2],it[0])}
    .collect().flatten().collate(3)

    // Get assembly data from snpdiffs
    snpdiff_assemblies = user_snpdiffs.map{it-> tuple(it[1],it[2])}
    .concat(user_snpdiffs.map{it-> tuple(it[10],it[11])})
    .join(exclude_ids,by:0,remainder:true)
    .filter{it -> it[0].toString() != "null"}
    .filter{it -> it[2].toString() != "Exclude"} 
    .map{it -> tuple(it[0],it[1],'SNPDiff')}
    .collect().flatten().collate(3)

    assembled_snpdiffs = snpdiff_assemblies
    .filter{it -> it[1].toString() != "null"}
    .unique{it->it[0]}.collect().flatten().collate(3)

    // Process any data provided as assemblies
    // Returns 2-item tuples with the following format: (Isolate_ID, Assembly_Path)
    ("${params.fasta}" != "" ? fetchQueryFasta() : Channel.empty()).set{query_fasta}
    ("${params.ref_fasta}" != "" ? fetchRefFasta() : Channel.empty()).set{ref_fasta}

    pre_assembled = assembled_snpdiffs
    .map{it -> tuple(it[0],it[1])}
    .concat(query_fasta)
    .concat(ref_fasta)
    .unique{it -> it[0]}
    .join(exclude_ids,by:0,remainder:true)
    .filter{it -> it[0].toString() != "null"}
    .filter{it -> it[2].toString() != "Exclude"}
    .map{it->tuple(it[0],it[1])}
    .collect().flatten().collate(2)

    // Process any data provided as reads
    // Returns 3-item tuples with the following format: (Isolate_ID, Read_Type, Read_Path)
    ("${params.reads}" != "" ? fetchQueryReads() : Channel.empty()).set{query_reads}
    ("${params.ref_reads}" != "" ? fetchRefReads() : Channel.empty()).set{ref_reads}

    all_reads = query_reads
    .concat(ref_reads)
    .unique{it->it[0]}
    .join(exclude_ids,by:0,remainder:true)
    .filter{it -> it[0].toString() != "null"}
    .filter{it -> it[3].toString() != "Exclude"}
    .map{it->tuple(it[0],it[1],it[2])}            
    .collect().flatten().collate(3)

    // Figure out if any assembly is necessary
    fasta_read_combo = all_reads.join(pre_assembled,by:0,remainder: true) |
    branch{it ->
        assembly: it[1].toString() == "null"
            return(tuple(it[0],it[2]))
        read: it[3].toString() == "null"
           return(tuple(it[0],it[1],it[2]))
        combo: true
           return(tuple(it[0],it[3]))}

    // Assemble reads if necessary
    assembled_reads = fasta_read_combo.read
    .collect().flatten().collate(3) | assembleReads

    // If runmode is 'assemble', tasks are complete
    if(params.runmode == "assemble"){
        query_data = Channel.empty()
        reference_data = Channel.empty()
    } else{
        
        // If FASTAs are provided via data and snpdiffs, use snpdiffs (as it's already been used)
        user_fastas = query_fasta
        .concat(ref_fasta)
        .concat(assembled_reads)
        .unique{it -> it[0]}
        .join(exclude_ids,by:0,remainder:true)
        .filter{it -> it[0].toString() != "null"}
        .filter{it -> it[2].toString() != "Exclude"}
        .map{it->tuple(it[0],it[1],'User')}
        .collect().flatten().collate(3)
        .join(assembled_snpdiffs,by:0,remainder:true)
        .filter{it -> it[3].toString() == "null"}
        .map{it->tuple(it[0],it[1])}

        // Get all assemblies
        all_assembled = assembled_snpdiffs
        .map{it -> tuple(it[0],it[1])}
        .concat(user_fastas)
        .unique{it->it[0]}.collect().flatten().collate(2)

        // Get data for isolates where a SNPDiff was provided, but no FASTA could be located
        no_assembly = snpdiff_assemblies
        .map{it -> tuple(it[0],it[1])}
        .filter{it -> it[1].toString() == "null"}
        .unique{it -> it[0]}
        .join(all_assembled,by:0,remainder:true)
        .filter{it -> it[2].toString() == "null"}
        .map{it->tuple(it[0],it[1])}
        .collect().flatten().collate(2)

        // Compile all samples
        all_samples = all_assembled
        .concat(no_assembly)
        .unique{it-> it[0]}.collect().flatten().collate(2)

        // If no reference data is provided return a blank channel
        if(!ref_mode){
            reference_data = Channel.empty()
            
            query_data = all_samples                
            .unique{it -> it[0]}
            .collect().flatten().collate(2)

        } else{

            // Process additional reference IDs
            ("${params.ref_id}" != "" ? processRefIDs() : Channel.empty()).set{user_ref_ids}
            
            all_ref_ids = ref_fasta.map{it->tuple(it[0])}
            .concat(ref_reads.map{it->tuple(it[0])})
            .concat(user_ref_ids)
            .unique{it-> it[0]}.collect().flatten().collate(1)
            .map{it -> tuple(it[0],"Reference")}
            .join(exclude_ids,by:0,remainder:true)
            .filter{it -> it[0].toString() != "null"}
            .filter{it -> it[2].toString() != "Exclude"}
            .map{it -> tuple(it[0],it[1])}

            reference_data = all_samples
            .join(all_ref_ids,by:0,remainder:true)
            .filter{it -> it[2].toString() == "Reference"}
            .map{it->tuple(it[0],it[1])}
            .unique{it -> it[0]}
            .collect().flatten().collate(2)

            // Save reference data to file
            reference_data
            .collect{it -> it[0]}
            | saveRefIDs

            if(params.runmode == "screen" || params.runmode == "align"){
                query_data = all_samples
                .join(all_ref_ids,by:0,remainder:true)
                .filter{it -> it[2].toString() != "Reference"}
                .map{it->tuple(it[0],it[1])}
                .unique{it -> it[0]}
                .collect().flatten().collate(2)
            } else if(params.runmode == "snp"){
                query_data = all_samples
                .unique{it -> it[0]}
                .collect().flatten().collate(2)
            }
        }
    }
}

// Fetching preassembled data //
workflow fetchQueryFasta{
    
    emit:
    query_fasta

    main:

    // If --fasta is set, grab assembly paths and characterize assemblies
    ("${params.fasta}" != "" ? getAssemblies(params.fasta) : Channel.empty()).set{query_fasta}
}
workflow fetchRefFasta{
    
    emit:
    ref_fasta

    main:

    // If --fasta is set, grab assembly paths and characterize assemblies
    ("${params.ref_fasta}" != "" ? getAssemblies(params.ref_fasta) : Channel.empty()).set{ref_fasta}
}
workflow getAssemblies{

    take:
    fasta_loc

    emit:
    fasta_data
    
    main:
    def trim_this = "${params.trim_name}"

    if(fasta_loc == ""){
        error "No assembly data provided via --fasta/--ref_fasta"
    } else{

        fasta_dir = file(fasta_loc)

        // If --fasta is a directory...
        if(fasta_dir.isDirectory()){
            ch_fasta = Channel.fromPath(["${fasta_dir}/*.fa","${fasta_dir}/*.fasta","${fasta_dir}/*.fna"])
        } 
        // If --fasta is a file...
        else if(fasta_dir.isFile()){
            
            // Check if it is a single fasta file...
            if(fasta_dir.getExtension() == "fa" || fasta_dir.getExtension() == "fna" || fasta_dir.getExtension() == "fasta"){
                ch_fasta = Channel.from(fasta_dir).map{it-> file(it)}
            } 
            // Otherwise, assume a file with paths to FASTAs
            else{
                ch_fasta = Channel.from(fasta_dir.readLines()).filter{ file -> file =~ /\.(fa|fasta|fna)$/}.map{it-> file(it)}
            }
        } else{
            error "$fasta_dir is not a valid directory or file..."
        }
        fasta_data = ch_fasta
        .filter { file(it).exists() }
        .map { filePath ->
            def fileName = file(filePath).getBaseName()
            def sampleName = fileName.replaceAll(trim_this, "")
            tuple(sampleName, filePath)}
    }
}
workflow processSNPDiffs{

    emit:
    snpdiffs_data
    
    main:

    if("${params.snpdiffs}" == ""){
        error "No assembly data provided via --snpdiffs"
    } else{

        snpdiffs_dir = file("${params.snpdiffs}")

        // If --fasta is a directory...
        if(snpdiffs_dir.isDirectory()){
            ch_snpdiffs = Channel.fromPath("${snpdiffs_dir}/*.snpdiffs")
        } 
        // If --fasta is a file...
        else if(snpdiffs_dir.isFile()){
            
            // Check if it is a single fasta file...
            if(snpdiffs_dir.getExtension() == "snpdiffs"){
                ch_snpdiffs = Channel.from(snpdiffs_dir)
            } 
            // Otherwise, assume a file with paths to SNPDiffs
            else{
                ch_snpdiffs = Channel.from(snpdiffs_dir.readLines()).filter{it->it.endsWith('.snpdiffs') }
            }
        } else{
            error "$snpdiffs_dir is not a valid directory or file..."
        }

        snpdiffs_data = ch_snpdiffs
            .filter { file(it).exists() }
            .collect() | getSNPDiffsData | splitCsv | collect | flatten | collate(19)

        // (1) SNPDiffs_File, (2) Query_ID, (3) Query_Assembly, (4) Query_Contig_Count, (5) Query_Assembly_Bases, 
        // (6) Query_N50, (7) Query_N90, (8) Query_L50, (9) Query_L90, (10) Query_SHA256,
        // (11) Reference_ID, (12) Reference_Assembly, (13) Reference_Contig_Count, (14) Reference_Assembly_Bases,
        // (15) Reference_N50, (16) Reference_N90, (17) Reference_L50, (18) Reference_L90, (19) Reference_SHA256
    }
}
process getSNPDiffsData{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(snpdiffs_paths)

    output:
    stdout

    script:

    user_snpdiffs_list.write(snpdiffs_paths.join('\n') + "\n")
    """
    $params.load_python_module
    python ${userSNPDiffs} --snpdiffs_file "${user_snpdiffs_list}" --trim_name "${params.trim_name}"
    """
}


// Fetching read data //
workflow fetchQueryReads{
    
    emit:
    query_reads

    main:

    // If --fasta is set, grab assembly paths and characterize assemblies
    ("${params.reads}" != "" ? processReads(params.reads,params.readext,params.forward,params.reverse) : Channel.empty()).set{query_reads}
}
workflow fetchRefReads{
    
    emit:
    ref_reads

    main:

    // If --fasta is set, grab assembly paths and characterize assemblies
    ("${params.ref_reads}" != "" ? processReads(params.ref_reads,params.ref_readext,params.ref_forward,params.ref_reverse) : Channel.empty()).set{ref_reads}
}
workflow processReads{

    take:
    read_loc
    read_ext
    forward
    reverse

    emit:
    read_info
    
    main:

    if(read_loc == ""){
        error "No data provided to --reads/--ref_reads"
    } else{

        read_dir = file(read_loc)

        // If --reads is a single directory, get all reads from that directory
        if(read_dir.isDirectory()){
            read_info = fetchReads(read_dir,read_ext,forward,reverse) | splitCsv
        } 

        // If --reads is a file including paths to many directories, process reads from all directories
        else if(read_dir.isFile()){
            read_info = fetchReads(Channel.from(read_dir.readLines()),read_ext,forward,reverse) | splitCsv
        }
        // Error if --reads doesn't point to a valid file or directory
        else{
            error "$read_dir is neither a valid file or directory..."
        }
    }      
}
process fetchReads{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val dir // Directory containing read files
    val read_ext // Extention for read files (e.g., fastq.gz or fq)
    val forward_suffix // Identifier for forward reads (e.g., _1.fastq or _R1_001.fq.gz)
    val reverse_suffix // Identifier for reverse reads (e.g., _2.fastq or _R2_001.fq.gz)

    output:
    stdout

    script:

    if(!file(dir).isDirectory()){
        error "$dir is not a valid directory..."
    } else{
    """
    $params.load_python_module
    python ${findReads} --read_dir ${dir} --read_filetype ${read_ext} --forward_suffix ${forward_suffix} --reverse_suffix ${reverse_suffix} --trim_name ${params.trim_name}
    """
    }
}

// Fetch reference IDs //
workflow processRefIDs{

    emit:
    ref_ids
    
    main:
    def trim_this = "${params.trim_name}"

    ref_ids = params.ref_id
    .tokenize(',')
    .unique()
    .collect { it ->
        "${it}".replaceAll(trim_this, "")}
    .flatten()
}

// Fetch reference IDs //
workflow processExclude{

    emit:
    exclude_ids
    
    main:
    def trim_this = "${params.trim_name}"

    exclude_ids = Channel.from(params.exclude
    .tokenize(',')
    .collect { it -> "${it}".replaceAll(trim_this, "")})
    .map{it -> tuple(it.toString(),"Exclude")}
    .unique{it -> it[0]}
}

process saveRefIDs{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(ref_ids)

    script:
    ref_id_file.append(ref_ids.join('\n') + '\n')        
    """
    """
}

// Assembly //
workflow assembleReads{

    take:
    to_assemble
    
    emit:
    assembled_data

    main:
    
    // Run SKESA on each entry
    assembly_output = skesaAssemble(to_assemble).splitCsv()

    // Print log of assemblies
    assembly_output.map {it -> it.join("\t")}.collect() | saveAssemblyLog

    // Return assembly data
    assembled_data = assembly_output.map{it->tuple(it[0],it[3])}
}
process skesaAssemble{
    memory '12 GB' // Add readcount/memory check?

    input:
    tuple val(sample_name),val(read_type),val(read_location)

    output:
    stdout
 
    script:
    assembly_file = file("${assembly_directory}/${sample_name}.fasta")
    
    // Ensure folder exists and file doesn't
    if(!assembly_directory.isDirectory()){
        error "$assembly_directory is not a valid directory..."
    } else if(assembly_file.isFile()){
        error "$assembly_file already exists..."
    } else if(read_type == "Paired"){
        forward_reverse = read_location.split(";")
        """
        $params.load_skesa_module
        skesa --cores ${skesa_cpus} --use_paired_ends --fastq ${forward_reverse[0]} ${forward_reverse[1]} --contigs_out ${assembly_file}
        echo "${sample_name},${read_type},${read_location},${assembly_file}"
        """
    } else if(read_type == "Single"){
        """
        $params.load_skesa_module
        skesa --cores ${skesa_cpus} --fastq ${read_location} --contigs_out ${assembly_file}
        echo "${sample_name},${read_type},${read_location},${assembly_file}"
        """
    } else{
        error "read_type should be Paired or Single, not $read_type..."
    }
}
process saveAssemblyLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(assembly_data)

    script:
    assembly_log.write(assembly_header)    
    assembly_log.append(assembly_data.join('\n') + '\n')    
    """
    """
}
