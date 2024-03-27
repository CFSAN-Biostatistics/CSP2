// Screening and SNP Pipeline processing
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)
screen_log_dir = file(params.screen_log_dir)
snp_log_dir = file(params.snp_log_dir)
screen_diffs_list = file("${log_directory}/Screen_Input.tsv")
screening_results_file = file(params.screening_results_file)

// Get QC thresholds
min_cov = params.min_cov.toFloat()
min_length = params.min_len.toInteger()
min_iden = params.min_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()

workflow runScreen {
    
    take:
    all_snpdiffs

    main:

    screening_results = all_snpdiffs
    .map{it -> it.join("\t")}
    .collect()
    | screenSNPDiffs
    | splitCsv() 

    screening_results
    .map{it -> it.join("\t")}
    .collect()
    | saveScreen
}

process screenSNPDiffs{

    input:
    val(all_snpdiffs)

    output:
    stdout

    script:

    screenDiffs = file("${projectDir}/bin/screenSNPDiffs.py")
    screen_diffs_list.write(all_snpdiffs.join('\n') + '\n')
    """
    $params.load_python_module
    $params.load_bedtools_module
    python $screenDiffs "${screen_diffs_list}" "${screen_log_dir}" "${min_cov}" "${min_length}" "${min_iden}" "${reference_edge}" "${query_edge}" "${params.dwin}" "${params.wsnps}" "${params.trim_name}"
    """
} 

process saveScreen{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(screen_results)

    script:
    screening_results_file.append(screen_results.join('\n') + '\n')
    """
    """
}


























workflow runSNPPipeline{
    take:
    snp_diff_data

    main:
    
    by_reference = snp_diff_data.map{tuple(it[1],it[2])}
    .groupTuple(by:0)
    .map { ref, diff_files -> tuple( ref.toString(), diff_files.collect() ) }
    | runSnpPipeline
}

process runSnpPipeline{

    input:
    tuple val(reference_id),val(diff_files)

    output:
    stdout

    script:

    out_dir = file(output_directory+"/SNP_${reference_id}")
    out_assembly = file(out_dir+"/SNPDiffs.txt")
    snp_script = file("${projectDir}/bin/runSNPPipeline.py")

    
    """
    $params.load_python_module
    $params.load_bedtools_module
    mkdir ${out_dir}
    cd ${out_dir}
    echo "${diff_files.join('\n')}" > $out_assembly
    python ${runSNP} "${reference_id}" "${out_assembly}" "${output_directory}"
    """
}


