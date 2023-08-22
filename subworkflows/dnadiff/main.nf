// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
params.out = "Yenta_${new java.util.Date().getTime()}"
params.outroot = ""
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

params.align_cov = 85
params.ref_iden = 99
params.ref_edge = 250
params.query_edge = 250
params.min_len = 500

alignment_coverage = params.align_cov.toFloat()
reference_identity = params.ref_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()
min_length = params.min_len.toInteger()

workflow runSnpPipeline{
    take:
    sample_data

    emit:
    sample_pairwise

    main:
    all_comparisons = sample_data.combine(sample_data).collect().flatten().collate(8).branch{
        
        same: "${it[0]}" == "${it[4]}"
        return(it)
        
        different: true
        return(it)}

    different_comparisons = all_comparisons.different.map{
        def lowerValue = "${it[0]}" <= "${it[4]}" ? "${it[0]}" : "${it[4]}"
        def higherValue = "${it[0]}" >= "${it[4]}" ? "${it[0]}" : "${it[4]}" 
        return tuple("${it[0]}","${it[4]}","${it[3]}","${it[7]}","${lowerValue};${higherValue}")}
        .collect().flatten().collate(5)
        .groupTuple(by:4).map{it -> tuple(it[2][0],it[2][1])}

    sample_pairwise = runMUmmer(different_comparisons) | splitCsv 
    | collect | flatten | collate(17)

    // Prep and save log files
    sample_log_file = prepSampleLog()
    sample_log_data = sample_data | join(sample_pairwise.map{it -> tuple(it[0],it[2],it[4])})
    saveSampleLog(sample_log_file,sample_log_data)

    snp_log_file = prepSNPLog()
    saveDNADiffLog(snp_log_file,sample_pairwise)

    // Run merging
    merged_snps = sample_pairwise.map { it.first() } | mergeSNPs
}

workflow runScreen{
    
    take:
    query_data
    reference_data

    main:
    
    // Make the pairwise comparisons of samples with references
    comparisons = query_data | combine(reference_data)
    | map{tuple(it[0],it[4],it[1],it[2],it[3],it[5],it[6],it[7])}
    // query,reference,query_datatype,query_data_location,query_fasta
    // reference_datatype,reference_data_location;reference_fasta

    // Run MUmmer jobs
    raw_screening_results = comparisons | map{tuple(it[4],it[7])} | runMUmmer | splitCsv
    // query,reference,query_seqs,ref_seqs,query_bases,
    // [5] ref_bases,percent_query_aligned_filtered,percent_ref_aligned_filtered,
    // [8] sample_category,final_snp_count,gsnps,rejected_snps_iden_count,
    // rejected_snps_edge_count,rejected_snps_dup_count,rejected_snps_density1000_count,
    // rejected_snps_density125_count,rejected_snps_density15_count

    // Prep log files
    query_log_file = prepQueryLog()
    ref_log_file = prepRefLog()
    dna_diff_log = prepDNADiffLog()

    // Save log files for isolates
    query_log_data = query_data | join(raw_screening_results.map{it -> tuple(it[0],it[2],it[4])})
    ref_log_data = reference_data | join(raw_screening_results.map{it -> tuple(it[1],it[3],it[5])})

    saveQueryLog(query_log_file,query_log_data)
    saveReferenceLog(ref_log_file,ref_log_data)
    
    // Save log files for DNADiff
    saveDNADiffLog(dna_diff_log,raw_screening_results)
}

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

process mergeSNPs{

    input:
    tuple val(query),val(reference), val(query_seqs), val(ref_seqs), val(query_bases), val(ref_bases), val(percent_query_aligned_filtered), val(percent_ref_aligned_filtered), val(sample_category), val(final_snp_count), val(gsnps), val(rejected_snps_iden_count), val(rejected_snps_edge_count), val(rejected_snps_dup_count), val(rejected_snps_density1000_count), val(rejected_snps_density125_count), val(rejected_snps_density15_count)

    output:
    val(snp_directory)

    script:
    """
    ${params.load_python_module}
    python $snp_script $mummer_directory $snp_directory
    """
}

// Log functions //
process prepQueryLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Query_Data.tsv"
    echo -n "${output_directory}/Query_Data.tsv"
    """
}
process prepSampleLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Isolate_Data.tsv"
    echo -n "${output_directory}/Isolate_Data.tsv"
    """
}
process prepRefLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    output:
    stdout

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Reference_Data.tsv"
    echo -n "${output_directory}/Reference_Data.tsv"
    """
}
process saveQueryLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(query_log_file)
    tuple val(sample_id),val(data_type),val(read_data),val(assembly_data),val(assembly_contigs),val(assembly_bases)

    script:

    if(!file(query_log_file).isFile()){
        error "$query_log_file doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${data_type}\t${read_data}\t${assembly_data}\t${assembly_contigs}\t${assembly_bases}" >> $query_log_file
    """
    }
}
process saveSampleLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(sample_log_file)
    tuple val(sample_id),val(data_type),val(read_data),val(assembly_data),val(assembly_contigs),val(assembly_bases)

    script:

    if(!file(sample_log_file).isFile()){
        error "$sample_log_file doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${data_type}\t${read_data}\t${assembly_data}\t${assembly_contigs}\t${assembly_bases}" >> $sample_log_file
    """
    }
}
process saveReferenceLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(reference_log_file)
    tuple val(sample_id),val(data_type),val(read_data),val(assembly_data),val(assembly_contigs),val(assembly_bases)

    script:

    if(!file(reference_log_file).isFile()){
        error "$reference_log_file doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${data_type}\t${read_data}\t${assembly_data}\t${assembly_contigs}\t${assembly_bases}" >> $reference_log_file
    """
    }
}
process prepDNADiffLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Query_ID\tReference_ID\tPercent_Reference_Covered\tPercent_Query_Covered\tYenta_SNPs\tCategory\tgSNPs\tFiltered_Edge\tFiltered_Identity\tFiltered_Duplicated\tRejected_Density_1000\tRejected_Density_125\tRejected_Density_15" > "${output_directory}/Screening_Results.tsv"
    echo -n "${output_directory}/Screening_Results.tsv"
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
process saveDNADiffLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(diff_log)
    tuple val(query),val(reference), val(query_seqs), val(ref_seqs), val(query_bases), val(ref_bases), val(percent_query_aligned_filtered), val(percent_ref_aligned_filtered), val(sample_category), val(final_snp_count), val(gsnps), val(rejected_snps_iden_count), val(rejected_snps_edge_count), val(rejected_snps_dup_count), val(rejected_snps_density1000_count), val(rejected_snps_density125_count), val(rejected_snps_density15_count)
    script:

    if(!file(diff_log).isFile()){
        error "$diff_log doesn't exist..."
    } else{
    """
    echo "${query}\t${reference}\t${percent_ref_aligned_filtered}\t${percent_query_aligned_filtered}\t${final_snp_count}\t${sample_category}\t${gsnps}\t${rejected_snps_edge_count}\t${rejected_snps_iden_count}\t${rejected_snps_dup_count}\t${rejected_snps_density1000_count}\t${rejected_snps_density125_count}\t${rejected_snps_density15_count}" >> $diff_log
    """
    }
}
