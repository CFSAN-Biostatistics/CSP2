// Subworkflow to fetch sample and reference data from --fasta/--reads/--ref_fasta/--ref_reads
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}
assembly_directory = file("${output_directory}/Assemblies")
assembly_log = file("${assembly_directory}/assembly.log")

// Set paths to accessory scripts
findPairedReads = file("${projectDir}/bin/fetchReads.py")
processFasta = file("${projectDir}/bin/processFasta.py")


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
                ch_fasta = Channel.from(fasta_dir)
            } 
            // Otherwise, assume a file with paths to FASTAs
            else{
                ch_fasta = Channel.from(fasta_dir.readLines())
            }
        } else{
            error "$fasta_dir is not a valid directory or file..."
        }
        fasta_data = ch_fasta 
        | processFasta 
        | splitCsv 
    }
}


// Fetching data that may require assembly //
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
            read_info = fetchPairedReads(read_dir,read_ext,forward,reverse) | splitCsv
        } 

        // If --reads is a file including paths to many directories, process reads from all directories
        else if(read_dir.isFile()){
            read_info = fetchPairedReads(Channel.from(read_dir.readLines()),read_ext,forward,reverse) | splitCsv
        }
        // Error if --reads doesn't point to a valid file or directory
        else{
            error "$read_dir is neither a valid file or directory..."
        }
    }      
}
process fetchPairedReads{

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
    python ${findPairedReads} ${dir} ${read_ext} ${forward_suffix} ${reverse_suffix} ${params.trim_name}
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

    // Create directory
    hold_file = prepAssemble() | collect | flatten
    
    // Run SKESA on each entry
    assembly_output = skesaAssemble(hold_file,to_assemble) | splitCsv

    // Print log of assemblies
    assembly_output | collect | flatten | collate(4) 
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}")} | collect 
    | saveAssemblyLog

    // Return assembly data
    assembled_data = assembly_output | map{it->it[3]} | processFasta | splitCsv
}
process prepAssemble{

    output:
    stdout

    script:
    if(assembly_directory.isDirectory()){
        error "$assembly_directory already exists..."
    } else{
        """
        mkdir $assembly_directory
        echo -n $assembly_directory
        """
    }
}
process skesaAssemble{

    // Set SKESA cores to 5 or fewer
    if("${params.cores}".toInteger() >= 5){
        cpus = 5
    } else{
        cpus = "${params.cores}".toInteger()
    }
    
    memory '6 GB'

    input:
    val (hold_file)
    tuple val(sample_name),val(read_type),val(read_location)

    output:
    stdout

    script:

    assembly_file = file("${assembly_directory}/${sample_name}.fasta")

    if(assembly_file.isFile()){
        error "$assembly_file already exists..."
    } else if(!file(assembly_directory).isDirectory()){
        error "$assembly_directory does not exist..."
    } else{
        if(read_type == "Paired"){
            forward_reverse = read_location.split(";")
            """
            $params.load_skesa_module
            skesa --use_paired_ends --fastq ${forward_reverse[0]} ${forward_reverse[1]} --contigs_out ${assembly_file}
            echo "$sample_name,$read_type,$read_location,$assembly_file"
            """
        } else if(read_type == "Single"){
            """
            $params.load_skesa_module
            skesa --fastq ${read_location} --contigs_out ${assembly_file}
            echo "$sample_name,$read_type,$read_location,$assembly_file"
            """
        } else{
            error "read_type should be Paired or Single, not $read_type..."
        }
    }
}
process saveAssemblyLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(assembly_data)

    script:
 
    """
    # Print header
    echo "Isolate_ID\tRead_Type\tRead_Data\tAssembly" > "${assembly_log}"
    echo "${assembly_data.join('\n')}" >> "${assembly_log}"
    """
}


// Process FASTA //
process processFasta{

    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val fasta_path

    output:
    stdout

    script:
    
    // Set path to accessory script
    fasta_script = file("${projectDir}/bin/processFasta.py")

    """
    module purge
    ${params.load_python_module}
    python ${fasta_script} ${fasta_path} ${params.trim_name}
    """
}




