// Subworkflow to fetch sample and reference data, assembling reads into genomes using SKESA if necessary
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}
assembly_directory = file("${output_directory}/Assemblies")

if(params.python_module == ""){
    params.load_python_module = ""
} else{
    params.load_python_module = "module load -s ${params.python_module}"
}
if(params.skesa_module == ""){
    params.load_skesa_module = ""
} else{
    params.load_skesa_module = "module load -s ${params.skesa_module}"
}

// Major workflows //
workflow fetchSampleData{
    
    emit:
    sample_data

    main:

    // Ensure data is provided for at least one sample
    if("${params.reads}" == "" && "${params.fasta}" == ""){
        error "No sample data specified by --reads/--fasta..."
    }

    // Fetch sample assemblies
    ("${params.fasta}" != "" ? getAssemblies(params.fasta) : Channel.empty()).set{sample_assembly_data}

    // Fetch sample reads
    ("${params.reads}" != "" ? getReads(params.reads,params.readext,params.forward,params.reverse) : Channel.empty()).set{sample_read_data}

    // Group by sample ID and identify samples where reads and assemblies are given
    grouped_data = sample_assembly_data.concat(sample_read_data) | collect | flatten | collate(3) | groupTuple |
    branch{it ->
        duo: it[1].size() == 2
            return(tuple(it[0], "Duo_"+it[1][1], it[2][1], it[2][0]))
        single: true
            if("${it[1][0]}" == "Assembly"){
                return tuple(it[0],it[1][0],"",it[2][0])
            } else{
                return tuple(it[0],it[1][0],it[2][0],"${assembly_directory}/${it[0]}.fasta")   
            }
    }

    assembled_data = grouped_data.single.filter { it -> it[1] == "Assembly" }
    unassembled_data = grouped_data.single.filter { it -> it[1] != "Assembly" } | skesaAssemble | splitCsv

    sample_data = grouped_data.duo.concat(assembled_data).concat(unassembled_data) | collect | flatten | collate(4)
}
workflow fetchReferenceData{

    take:
    ref_reads
    ref_fasta
    
    emit:
    reference_data

    main:
        
    if(ref_reads == "" && ref_fasta == ""){
        reference_data = Channel.empty()
    } else{
        // Collect paths to read/assembly data for references
        ("${ref_reads}" != "" ? getReads(ref_reads,params.ref_readext,params.ref_forward,params.ref_reverse) : Channel.empty()).set{reference_read_data}
        ("${ref_fasta}" != "" ? getAssemblies(ref_fasta) : Channel.empty()).set{reference_assembly_data}

        // Group by sample ID and identify samples where reads and assemblies are given
        grouped_data = reference_assembly_data.concat(reference_read_data) | collect | flatten | collate(3) | groupTuple |
        branch{it ->
            duo: it[1].size() == 2
                return(tuple(it[0], "Duo_"+it[1][1], it[2][1], it[2][0]))
            single: true
                if("${it[1][0]}" == "Assembly"){
                    return tuple(it[0],it[1][0],"",it[2][0])
                } else{
                    return tuple(it[0],it[1][0],it[2][0],"${assembly_directory}/${it[0]}.fasta")   
                }
        }

        assembled_data = grouped_data.single.filter { it -> it[1] == "Assembly" }
        unassembled_data = grouped_data.single.filter { it -> it[1] != "Assembly" } | skesaAssemble | splitCsv

        reference_data = grouped_data.duo.concat(assembled_data).concat(unassembled_data) | collect | flatten | collate(4)
    }
}

// Fetch workflows //
workflow getReads{

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
            read_info = fetchPairedReads(read_dir,read_ext,forward,reverse) 
            | splitCsv 
            | map{tuple(it[0].toString(),it[1].toString(),it[2].toString())}
        } 

        // If --reads is a file including paths to many directories, process reads from all directories
        else if(read_dir.isFile()){
            read_info = fetchPairedReads(Channel.from(read_dir.readLines()),read_ext,forward,reverse) 
            | splitCsv 
            | map{tuple(it[0].toString(),it[1].toString(),it[2].toString())}

        }
        // Error if --reads doesn't point to a valid file or directory
        else{
            error "$read_dir is neither a valid file or directory..."
        }
    }      
}
workflow getAssemblies{

    take:
    fasta_loc

    emit:
    fasta_data
    
    main:

    if(fasta_loc == ""){
        error "No assembly data provided via --fasta"
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
        }
        else{
            error "$fasta_dir is not a valid directory or file..."
        }

        fasta_data = ch_fasta
        | map{tuple(file("$it").getBaseName(),"Assembly",file("$it"))} // Get sample name from filename, return tuple of ID, 'Assembly',fasta location
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
    
    // Set path to accessory script
    findPairedReads = file("${projectDir}/bin/fetchReads.py")

    """
    ${params.load_python_module}
    python ${findPairedReads} ${dir} ${read_ext} ${forward_suffix} ${reverse_suffix}
    """
}

// Assembly //

process skesaAssemble{

    // Set SKESA cores to 5 or fewer
    if("${params.cores}".toInteger() >= 5){
        cpus = 5
    } else{
        cpus = "${params.cores}".toInteger()
    }
    
    memory '6 GB'

    input:
    tuple val(sample_name),val(read_type),val(read_location),val(assembly_out)

    output:
    stdout

    script:

    assembly_file = file(assembly_out)
    assembly_dir = assembly_file.getParent()

    if(assembly_file.isFile()){
        error "$assembly_out already exists..."
    } else if(!assembly_dir.isDirectory()){
        error "$assembly_dir does not exist..."
    } else{
        if(read_type == "Paired"){
            forward_reverse = read_location.split(";")
            """
            $params.load_skesa_module
            skesa --use_paired_ends --fastq ${forward_reverse[0]} ${forward_reverse[1]} --contigs_out ${assembly_file}
            echo "$sample_name,$read_type,$read_location,$assembly_out"
            """
        } else if(read_type == "Single"){
            """
            $params.load_skesa_module
            skesa --fastq ${read_location} --contigs_out ${assembly_file}
            echo "$sample_name,$read_type,$read_location,$assembly_out"
            """
        } else{
            error "read_type should be Paired or Single, not $read_type..."
        }
    }
}
