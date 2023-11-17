// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

log_directory = file("${output_directory}/logs")
snp_diffs_file = file("${log_directory}/SNPDiffs.txt")

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

// Set paths to accessory scripts
screenDiffs = file("${projectDir}/bin/screenSNPDiffs.py")
runSNP = file("${projectDir}/bin/runSNPPipeline.py")

// Set modules
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
params.load_bedtools_module = params.bedtools_module == "" ? "" : "module load -s ${params.bedtools_module}"


workflow runScreen{
    take:
    snp_diff_data

    main:

    snp_diff_data.collect{it[2]} | screenSNPDiffs
}

workflow runSNPPipeline{
    take:
    snp_diff_data

    main:
    
    by_reference = snp_diff_data.map{tuple(it[1],it[2])}
    .groupTuple(by:0)
    .map { ref, diff_files -> tuple( ref.toString(), diff_files.collect() ) }
    | runSnpPipeline

    by_reference.subscribe{println("$it")}


}

process runSnpPipeline{

    input:
    tuple val(reference_id),val(diff_files)

    output:
    stdout

    script:

    out_dir = file(output_directory+"/SNP_${reference_id}")
    out_assembly = file(out_dir+"/SNPDiffs.txt")
    
    """
    $params.load_python_module
    $params.load_bedtools_module
    mkdir ${out_dir}
    cd ${out_dir}
    echo "${diff_files.join('\n')}" > $out_assembly
    python ${runSNP} "${reference_id}" "${out_assembly}" "${output_directory}"
    """
}

process screenSNPDiffs{

    input:
    val(snp_diffs)

    script:
    
    """
    $params.load_python_module
    echo "${snp_diffs.join('\n')}" > $snp_diffs_file
    python ${screenDiffs} "${snp_diffs_file}" "${output_directory}"
    """
}

