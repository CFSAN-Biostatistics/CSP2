// Subworkflow to fetch sample and reference data from --fasta/--reads/--ref_fasta/--ref_reads

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

// Save assembly data in the main directory if --runmode is 'assemble'
if(run_mode == "assemble"){
    assembly_directory = file("${output_directory}")
    log_directory = file("${output_directory}")
} else{
    log_directory = file("${output_directory}/logs")
    assembly_directory = file("${output_directory}/Assemblies")
}
assembly_log = file("${log_directory}/Assembly_Data.tsv")
user_snpdiffs_list = file("${log_directory}/Imported_SNPDiffs.txt")

// Set paths to accessory scripts
findReads = file("${projectDir}/bin/fetchReads.py")
processFasta = file("${projectDir}/bin/processFasta.py")
userSNPDiffs = file("${projectDir}/bin/userSNPDiffs.py")

// Set up modules if needed
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
params.load_skesa_module = params.skesa_module == "" ? "" : "module load -s ${params.skesa_module}"

// Set SKESA cores to 5 or fewer
skesa_cpus = (params.cores as Integer) >= 5 ? 5 : (params.cores as Integer)

// Top-level workflow //
workflow fetchData{

    emit:
    query_data
    reference_data
    snpdiff_data

    main:

    // Process snpdiffs alignments
    ("${params.snpdiffs}" != "" ? processSNPDiffs() : Channel.empty()).set{user_snpdiffs}
    snpdiff_data = user_snpdiffs.map{it -> tuple(it[0],it[2],it[4])}.collect().flatten().collate(3)

    // Grab assemblies from snpdiffs if they still exist in the same directory
    snpdiff_assemblies = user_snpdiffs.map{it-> tuple(it[0],it[1])}
    .concat(user_snpdiffs.map{it-> tuple(it[2],it[3])})
    .filter{it -> it[1] != null && it[1].toString() != "null"}
    .unique().collect().flatten().collate(2)

    // Process any data provided as assemblies
    ("${params.fasta}" != "" ? fetchQueryFasta() : Channel.empty()).set{query_fasta}
    ("${params.ref_fasta}" != "" ? fetchRefFasta() : Channel.empty()).set{ref_fasta}

    all_assembled = snpdiff_assemblies.concat(query_fasta).concat(ref_fasta).unique().collect().flatten().collate(2)

    // Process any data provided as reads
    ("${params.reads}" != "" ? fetchQueryReads() : Channel.empty()).set{query_reads}
    ("${params.ref_reads}" != "" ? fetchRefReads() : Channel.empty()).set{ref_reads}
    all_reads =  query_reads.concat(ref_reads).collect().flatten().collate(3).unique{it->it[0]}
    
    // Figure out if any assembly is necessary
    fasta_read_combo = all_reads.join(all_assembled,by:0,remainder: true) |
    branch{it ->
        assembly: it[1].toString() == "null"
            return(tuple(it[0],it[2]))
        read: it[3].toString() == "null"
           return(tuple(it[0],it[1],it[2]))
        combo: true
           return(tuple(it[0],it[3]))}

    assembled_reads = fasta_read_combo.read.collect().flatten().collate(3).unique{it->it[0]} | assembleReads
    user_fastas = query_fasta.concat(ref_fasta).concat(assembled_reads).unique{it->it[0]}.collect().flatten().collate(2)

    assembled_isolates = snpdiff_assemblies.map{it-> tuple(it[0],it[1])}
    .concat(user_fastas.map{it -> tuple(it[0],it[1],'User')}.join(snpdiff_assemblies.map{it-> tuple(it[0],it[1],"SNPDiff")},by:0,remainder:true)
    .filter( (it->it[3] == null || it[3].toString() == "null"))
    .map{it-> tuple(it[0],it[1])}).collect().flatten().collate(2)

    // Process additional reference IDs
    ("${params.ref_id}" != "" ? processRefIDs() : Channel.empty()).set{user_ref_ids}

    reference_data = ref_fasta.map{it->tuple(it[0])}
    .concat(ref_reads.map{it->tuple(it[0])})
    .concat(user_ref_ids)
    .unique{it-> it[0]}.collect().flatten().collate(2)
    .join(assembled_isolates,by:0).collect().flatten().collate(2)

    if((run_mode == "assemble") || (run_mode == "snp")){
        query_data = assembled_isolates
    } else if ((run_mode == "align") || (run_mode == "screen")){
        query_data = assembled_isolates.map{it->tuple(it[0],it[1],"Query")}
        .join(reference_data.map{it->tuple(it[0],it[1],"Reference")},by:0,remainder:true)
        .filter( (it->it[3] == null || it[3].toString() == "null"))
        .map{it-> tuple(it[0],it[1])}
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

    def trim_this = "${params.trim_name}"

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
            .collect() | getSNPDiffsData | splitCsv | collect | flatten | collate(5)
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
    """
    ${params.load_python_module}

    echo "${snpdiffs_paths.join('\n')}" > $user_snpdiffs_list
    python ${userSNPDiffs} "${user_snpdiffs_list}" ${params.trim_name}
    """
}

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
    ${params.load_python_module}
    python ${findReads} ${dir} ${read_ext} ${forward_suffix} ${reverse_suffix} ${params.trim_name}
    """
    }
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
        $params.load_python_module
        $params.load_skesa_module
        skesa --cores ${skesa_cpus} --use_paired_ends --fastq ${forward_reverse[0]} ${forward_reverse[1]} --contigs_out ${assembly_file}
        python ${processFasta} "${sample_name}" "${read_type}" "${read_location}" "${assembly_file}"
        """
    } else if(read_type == "Single"){
        """
        $params.load_python_module
        $params.load_skesa_module
        skesa --cores ${skesa_cpus} --fastq ${read_location} --contigs_out ${assembly_file}
        python ${processFasta} "${sample_name}" "${read_type}" "${read_location}" "${assembly_file}"
        """
    } else{
        error "read_type should be Paired or Single, not $read_type..."
    }
}

// Logging Processes//
process saveAssemblyLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(assembly_data)

    script:
 
    """
    echo "Isolate_ID\tRead_Type\tRead_Data\tAssembly\tContig_Count\tAssembly_Bases\tN50\tL50\tN90\tL90\tSHA256" > "${assembly_log}"
    echo "${assembly_data.join('\n')}" >> "${assembly_log}"
    """
}