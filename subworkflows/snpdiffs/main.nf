// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['all','assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'all', 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

log_directory = file("${output_directory}/logs")
mummer_directory = file("${output_directory}/MUMmer_Output")
snpdiffs_directory = file("${output_directory}/snpdiffs")
snpdiffs_list_file = file("${log_directory}/All_SNPDiffs.txt")
ref_id_file = file("${log_directory}/Ref_IDs.txt")

snp_directory = file("${output_directory}/SNP_Analysis")

// Set paths to accessory scripts
screen_script = file("${projectDir}/bin/screenSNPDiffs.py")
snp_script = file("${projectDir}/bin/runSNPPipeline.py")

// Set modules
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
params.load_bedtools_module = params.bedtools_module == "" ? "" : "module load -s ${params.bedtools_module}"

// Get QC thresholds
min_cov = params.min_cov.toFloat()
min_length = params.min_len.toInteger()
min_iden = params.min_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()

workflow runScreen {
    
    take:
    all_snpdiffs
    reference_data

    main:

    ref_ids = reference_data.collect{it[0]}
    snpdiff_files = all_snpdiffs.collect{it[2]}
    screenSNPDiffs(snpdiff_files,ref_ids)
}

process screenSNPDiffs{

    input:
    val(snp_diffs)
    val(ref_ids)

    script:
    if((params.ref_fasta == "") && (params.ref_reads == "") && (params.ref_id == "")){
    """
    $params.load_python_module
    echo "${snp_diffs.join('\n')}" > $snp_diffs_file
    python ${screenDiffs} "${snp_diffs_file}" "${output_directory}"
    """
    } else{
    """
    $params.load_python_module
    echo "${snp_diffs.join('\n')}" > $snp_diffs_file
    echo "${ref_ids.join('\n')}" > $ref_id_file
    python ${screenDiffs} "${snp_diffs_file}" "${output_directory}"
    """   
    }
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
    
    """
    $params.load_python_module
    $params.load_bedtools_module
    mkdir ${out_dir}
    cd ${out_dir}
    echo "${diff_files.join('\n')}" > $out_assembly
    python ${runSNP} "${reference_id}" "${out_assembly}" "${output_directory}"
    """
}


