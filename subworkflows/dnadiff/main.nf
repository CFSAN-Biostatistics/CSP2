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

///// SNP PIPELINE //////
workflow runSnpPipeline{
    take:
    sample_data
    reference_data

    main:
    comparisons = sample_data.combine(reference_data).collect().flatten().collate(8)
    .filter{it[0] != it[4]}
    .map{
        def lowerValue = "${it[0]}" < "${it[4]}" ? "${it[0]}" : "${it[4]}"
        def higherValue = "${it[0]}" > "${it[4]}" ? "${it[0]}" : "${it[4]}" 
        return tuple("${it[0]}","${it[4]}","${it[3]}","${it[7]}","${lowerValue};${higherValue}")}
        .collect().flatten().collate(5).map{it -> tuple(it[2],it[3])}
    
    sample_pairwise = runMUmmer(comparisons) | splitCsv 
    | collect | flatten | collate(18)

    // Save sample log
    log_data_a = sample_data | join(sample_pairwise.map{it -> tuple(it[0],it[2],it[4])})
    log_data_b = sample_data | join(sample_pairwise.map{it -> tuple(it[1],it[3],it[5])})
    
    log_data_a.concat(log_data_b) | toSortedList({ a, b -> a[0] <=> b[0] }) | collect | flatten | collate(6) | distinct
    | map { it -> tuple(it[0].toString(),it[1].toString(),it[2].toString(),it[3].toString(),it[4].toString(),it[5].toString())}
    | collect | flatten | collate(6)
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\t${it[4]}\t${it[5]}")}
    | collect
    | saveSampleLog

    // Save DNA diff results, return true when done
    diff_results = sample_pairwise
    | map { it -> tuple(it[0].toString(),it[1].toString(),it[6].toString(),it[7].toString(),it[8].toString(),it[9].toString(),it[10].toString(),it[11].toString(),it[12].toString(),it[13].toString(),it[14].toString(),it[15].toString(),it[16].toString(),it[17].toString())}
    | collect | flatten | collate(14)
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\t${it[4]}\t${it[5]}\t${it[6]}\t${it[7]}\t${it[8]}\t${it[9]}\t${it[10]}\t${it[11]}\t${it[12]}\t${it[13]}")}
    | collect
    | saveDNADiffLog | collect | flatten | first 

    // Run merging + tree building after log is written
    merged_snps = refSNPs(diff_results,reference_data)
}
process refSNPs{

    input:
    val(ready)
    tuple val(ref_isolate),val(sample_type),val(read_location),val(assembly_location)
    
    output:
    val(snp_directory)

    script:
    """
    ${params.load_python_module}
    ${params.load_bedtools_module}
    python $ref_snp_script $output_directory $alignment_coverage $reference_identity $min_length $ref_isolate
    """
}

////// SCREENER //////
workflow runScreen{
    
    take:
    query_data
    reference_data

    main:
    
    // Make the pairwise comparisons of samples with references
    comparisons = query_data | combine(reference_data)
    | map{tuple(it[0],it[4],it[1],it[2],it[3],it[5],it[6],it[7])}
    // query,reference,
    // query_datatype,query_data_location,query_fasta
    // reference_datatype,reference_data_location;reference_fasta

    // Run MUmmer jobs
    raw_screening_results = comparisons | map{tuple(it[4],it[7])} | runMUmmer | splitCsv | collect | flatten | collate(18)
    // query,reference,query_seqs,ref_seqs,query_bases,
    // [5] ref_bases,percent_query_aligned_filtered,percent_ref_aligned_filtered,
    // [8] sample_category,final_snp_count,gsnps,rejected_snps_iden_count,
    // rejected_snps_edge_count,rejected_snps_dup_count,rejected_snps_density1000_count,
    // rejected_snps_density125_count,rejected_snps_density15_count

    // Save log files for queries
    query_data | join(raw_screening_results.map{it -> tuple(it[0],it[2],it[4])}) 
    | map { it -> tuple(it[0].toString(),it[1].toString(),it[2].toString(),it[3].toString(),it[4].toString(),it[5].toString())}
    | collect | flatten | collate(6)
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\t${it[4]}\t${it[5]}")}
    | collect
    | saveQueryLog

    // Save log files for references    
    ref_log_data = reference_data | join(raw_screening_results.map{it -> tuple(it[1],it[3],it[5])})
    | map { it -> tuple(it[0].toString(),it[1].toString(),it[2].toString(),it[3].toString(),it[4].toString(),it[5].toString())}
    | collect | flatten | collate(6)
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\t${it[4]}\t${it[5]}")}
    | collect
    | saveReferenceLog
    
    // Save log files for DNADiff
    raw_screening_results | map { it -> tuple(it[0].toString(),it[1].toString(),it[6].toString(),it[7].toString(),it[8].toString(),it[9].toString(),it[10].toString(),it[11].toString(),it[12].toString(),it[13].toString(),it[14].toString(),it[15].toString(),it[16].toString(),it[17].toString())}
    | collect | flatten | collate(14)
    | map { it -> tuple("${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\t${it[4]}\t${it[5]}\t${it[6]}\t${it[7]}\t${it[8]}\t${it[9]}\t${it[10]}\t${it[11]}\t${it[12]}\t${it[13]}")}
    | collect
    | saveScreeningLog
}

///// MUMmer //////
process runMUmmer{

    newForks = "${params.cores}".toInteger() * "${params.maxForks}".toInteger()
    maxForks = newForks.toInteger()
    cpus = 1
    
    input:
    tuple val(query_fasta),val(ref_fasta)
    
    output:
    stdout

    script:

    query_name = file(query_fasta).getBaseName()
    ref_name = file(ref_fasta).getBaseName()
    report_id = "${query_name}_vs_${ref_name}"

    // Ensure MUmmer directories exist
    if(!mummer_directory.isDirectory()){
        error "$mummer_directory does not exist..."
    } else if(!raw_mummer_directory.isDirectory()){
        error "$raw_mummer_directory does not exist..."
    } else{
        """
        $params.load_mummer_module
        $params.load_python_module
        $params.load_bedtools_module

        cd ${raw_mummer_directory}
        dnadiff -p ${report_id} ${ref_fasta} ${query_fasta}
        python ${mummer_processing_script} ${query_name} ${ref_name} ${report_id} ${raw_mummer_directory} ${alignment_coverage} ${reference_identity} ${reference_edge} ${query_edge} ${min_length}       
        rm -rf ${raw_mummer_directory}/${report_id}.mdelta
        rm -rf ${raw_mummer_directory}/${report_id}.mcoords
        rm -rf ${raw_mummer_directory}/${report_id}.1delta
        rm -rf ${raw_mummer_directory}/${report_id}.delta
        rm -rf ${raw_mummer_directory}/${report_id}.qdiff
        rm -rf ${raw_mummer_directory}/${report_id}.rdiff
        rm -rf ${raw_mummer_directory}/${report_id}.unref
        rm -rf ${raw_mummer_directory}/${report_id}.unqry
        mv ${raw_mummer_directory}/${report_id}.snps ${raw_mummer_directory}/SNPs
        mv ${raw_mummer_directory}/${report_id}.report ${raw_mummer_directory}/Reports
        mv ${raw_mummer_directory}/${report_id}.1coords ${raw_mummer_directory}/1coords
        """
    }
}

// Log functions //
process saveSampleLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(log_data)

    script:
 
    """
    # Print header
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Isolate_Data.tsv"
    echo "${log_data.join('\n')}" >> "${output_directory}/Isolate_Data.tsv"
    """
}
process saveDNADiffLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(diff_data)

    output:
    val(true)

    script:

    """
    echo "Query_ID\tReference_ID\tPercent_Reference_Covered\tPercent_Query_Covered\tCategory\tYenta_SNPs\tMedian_SNP_Perc_Iden\tgSNPs\tFiltered_Edge\tFiltered_Identity\tFiltered_Duplicated\tRejected_Density_1000\tRejected_Density_125\tRejected_Density_15" > "${output_directory}/Raw_Pairwise_Distances.tsv"
    echo "${diff_data.join('\n')}" >> "${output_directory}/Raw_Pairwise_Distances.tsv"
    """
}
process saveScreeningLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(diff_data)

    script:

    """
    echo "Query_ID\tReference_ID\tPercent_Reference_Covered\tPercent_Query_Covered\tCategory\tYenta_SNPs\tMedian_SNP_Perc_Iden\tgSNPs\tFiltered_Edge\tFiltered_Identity\tFiltered_Duplicated\tRejected_Density_1000\tRejected_Density_125\tRejected_Density_15" > "${output_directory}/Raw_Pairwise_Distances.tsv"
    echo "${diff_data.join('\n')}" >> "${output_directory}/Reference_Screening_Results.tsv"
    """
}
process saveQueryLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(log_data)

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Query_Data.tsv"
    echo "${log_data.join('\n')}" >> "${output_directory}/Query_Data.tsv"
    """
}
process saveReferenceLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(log_data)

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Reference_Data.tsv"
    echo "${log_data.join('\n')}" >> "${output_directory}/Reference_Data.tsv"
    """
}
