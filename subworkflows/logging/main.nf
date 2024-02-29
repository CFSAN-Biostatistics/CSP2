// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

// Assess run mode
if (params.runmode == "") {
    if (params.snpdiffs != "") {
        run_mode = "screen" // If .snpdiffs are provided, generate a summary of query/reference alignments
    } else if((params.reads != "" || params.ref_reads != "") && (params.fasta == "" && params.ref_fasta == "")){
        run_mode = "assemble" // If only reads are provided, generate assemblies
    } else if(params.reads == "" && params.fasta == ""){
        error "No query data provided via --reads/--fasta/--snpdiffs" // Exit if no data is provided
    } else if (params.ref_fasta == "" && params.ref_reads == ""){
        run_mode = "snp" // If query data is provided without reference data, run the SNP pipeline with RefChooser
    } else if((params.reads != "" || params.fasta != "") && (params.ref_reads != "" || params.ref_fasta == "")){
        run_mode = "screen" // If query and reference data are provided, perform MUMmer alignment and generate a summary
    } else if((params.fasta == "" && params.reads == "") && (params.ref_fasta != "" || params.ref_reads != "")){
        error "Reference data provided via --ref_reads/--ref_fasta, but no query data provided by --reads/--fasta/--snpdiffs" // Exit if no query data is provided
    } 
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Save assembly data in the main directory if --runmode is 'assemble'
if(run_mode == "assemble"){
    assembly_directory = file("${output_directory}")
    assembly_log = file("${output_directory}/Isolate_Data.tsv")
} else{
    log_directory = file("${output_directory}/logs")
    assembly_directory = file("${output_directory}/Assemblies")
    assembly_log = file("${log_directory}/Assembly_Log.tsv")
    isolate_file = file("${output_directory}/Isolate_Data.tsv")
}

// Logging Processes//
process saveIsolateLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(isolate_data)

    output:
    val(isolate_data)

    script:
    """
    echo "Isolate_ID\tRead_Type\tRead_Data\tAssembly\tContig_Count\tAssembly_Bases\tN50\tL50\tN90\tL90\tSHA256" > "${isolate_file}"
    echo "${isolate_data.join('\n')}" >> "${isolate_file}"
    """
}

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