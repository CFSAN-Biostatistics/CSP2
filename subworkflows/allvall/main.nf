// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}

mummer_directory = file("${output_directory}/MUmmer_Output")
snp_directory = file("${output_directory}/SNP_Analysis")
raw_mummer_directory = file("${output_directory}/MUmmer_Output/Raw")

// Set path to mummer script
mummer_processing_script = file("$projectDir/bin/filterMUmmer.py")

// Set path to snp script
snp_script = file("$projectDir/bin/mergeSNPs.py")
ref_snp_script = file("$projectDir/bin/refSNPs.py")

// Set modules if necessary
if(params.mummer_module == ""){
    params.load_mummer_module = ""
} else{
    params.load_mummer_module = "module load -s ${params.mummer_module}"
}
params.python_module = ""
if(params.python_module == ""){
    params.load_python_module = ""
} else{
    params.load_python_module = "module load -s ${params.python_module}"
}
params.bedtools_module = ""
if(params.bedtools_module == ""){
    params.load_bedtools_module = ""
} else{
    params.load_bedtools_module = "module load -s ${params.bedtools_module}"
}

alignment_coverage = params.align_cov.toFloat()
reference_identity = params.ref_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()
min_length = params.min_len.toInteger()

workflow runAllvAll{
    take:
    sample_data

    main:

    comparisons = sample_data.combine(sample_data).collect().flatten().collate(8)
    .filter{it[0] != it[4]}
    .map{
        def lowerValue = "${it[0]}" < "${it[4]}" ? "${it[0]}" : "${it[4]}"
        def higherValue = "${it[0]}" > "${it[4]}" ? "${it[0]}" : "${it[4]}" 
        return tuple("${it[0]}","${it[4]}","${it[3]}","${it[7]}","${lowerValue};${higherValue}")}
    .collect().flatten().collate(5)
    .groupTuple(by:4).map{it -> tuple(it[2][0],it[2][1])}

    sample_pairwise = runMUmmer(comparisons) | splitCsv 
    | collect | flatten | collate(17)

    // Prep and save log files
    sample_log_file = prepSampleLog()
    log_data_a = sample_data | join(sample_pairwise.map{it -> tuple(it[0],it[2],it[4])})
    log_data_b = sample_data | join(sample_pairwise.map{it -> tuple(it[1],it[3],it[5])})
    sample_log_data = log_data_a.concat(log_data_b) | toSortedList({ a, b -> a[0] <=> b[0] }) | flatten | collate(6) | distinct
    saveSampleLog(sample_log_file,sample_log_data)

    // Move this to after merging?
    snp_log_file = prepSNPLog()
    diff_results = saveDNADiffLog(snp_log_file,sample_pairwise)

    // Run merging + tree building
    merged_snps = diff_results | collect | mergeSNPs
}

process mergeSNPs{

    input:
    val(ready)
    
    output:
    val(snp_directory)

    script:
    """
    ${params.load_python_module}
    python $snp_script $output_directory $mummer_directory $snp_directory $alignment_coverage $perc_max_n
    """
}

process prepSNPLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Query_ID\tReference_ID\tPercent_Reference_Covered\tPercent_Query_Covered\tYenta_SNPs\tCategory\tgSNPs\tFiltered_Edge\tFiltered_Identity\tFiltered_Duplicated\tRejected_Density_1000\tRejected_Density_125\tRejected_Density_15" > "${output_directory}/Raw_Pairwise_Distances.tsv"
    echo -n "${output_directory}/Raw_Pairwise_Distances.tsv"
    """
}