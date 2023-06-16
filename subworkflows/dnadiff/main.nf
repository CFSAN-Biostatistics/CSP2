// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output folder
params.out = "./YENTA_${new java.util.Date().getTime()}"
output_directory = file("${params.out}")
results_directory = file("${output_directory}/Screening_Results")
reference_directory = file("${output_directory}/Reference_Strain_Data")

// Set MUmmer module if necessary
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
params.query_edge = 500

alignment_coverage = params.align_cov.toFloat()
reference_identity = params.ref_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()

workflow dnaDiff{
    
    take:
    sample_data
    reference_data

    main:

    // Prep log files
    getSampleLog = prepSampleLog()
    getReferenceLog = prepReferenceLog()
    getDNADiffLog = prepDNADiffLog()

    comparisons = sample_data
    | combine(reference_data) // Make the pairwise comparisons of samples with references
    | runMUmmer // Run MUmmer jobs
    | splitCsv

    // Wait for summaries to finish before making logs
    log_data = comparisons | collect | flatten | collate(23)
    
    // Save log files for samples
    sample_log_data = log_data.map { it[0] } | flatten | unique | join(log_data.map{it -> tuple(it[0],it[1],it[2],it[3],it[4],it[5])})
    saveSampleLog(getSampleLog,sample_log_data)

    // Save log files for reference
    reference_log_data = log_data.map { it[6] } | flatten | unique | join(log_data.map{it -> tuple(it[6],it[7],it[8],it[9],it[10],it[11])})
    saveReferenceLog(getReferenceLog,reference_log_data)
    
    // Save log files for DNADiff
    saveDNADiffLog(getDNADiffLog,comparisons.map{tuple(it[0],it[6],it[12],it[13],it[14],it[15],it[16],it[17],it[18],it[19],it[20],it[21],it[22])})
}

process runMUmmer{
     
    input:
    tuple val(isolate),val(data_type),val(data_location),val(isolate_fasta),val(ref_strain),val(ref_data_type),val(ref_data_location),val(ref_fasta)
    
    output:
    stdout

    script:

    mummer_dir = file("${output_directory}/${isolate}/MUmmer")
    mummer_filter_script = "$projectDir/bin/filterMUmmer.py"

    // Replace blank with NA if assemblies are provided
    if(data_location == ""){
        data_location = "NA"
    }
    if(ref_data_location == ""){
        ref_data_location = "NA"
    }
    if(!mummer_dir.isDirectory()){
        error "$mummer_dir does not exist..."
    } else{
        """
        $params.load_mummer_module
        $params.load_python_module
        $params.load_bedtools_module

        cd $mummer_dir  
        dnadiff -p ${ref_strain} ${ref_fasta} ${isolate_fasta}
        python ${mummer_filter_script} ${isolate} ${data_type} "${data_location}" ${isolate_fasta} ${ref_strain} ${ref_data_type} "${ref_data_location}" ${ref_fasta} ${mummer_dir} ${alignment_coverage} ${reference_identity} ${reference_edge} ${query_edge}        
        """
    }
}

// Log functions //
process prepSampleLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Sample_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > ${output_directory}/Sample_Data.tsv
    echo -n "${output_directory}/Sample_Data.tsv"
    """
}
process prepReferenceLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    output:
    stdout

    script:
    """
    echo "Reference_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > ${reference_directory}/Reference_Data.tsv
    echo -n "${reference_directory}/Reference_Data.tsv"
    """
}
process saveSampleLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(sample_log)
    tuple val(sample_id),val(data_type),val(read_data),val(assembly_data),val(assembly_contigs),val(assembly_bases)

    script:

    if(!file(sample_log).isFile()){
        error "$sample_log doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${data_type}\t${read_data}\t${assembly_data}\t${assembly_contigs}\t${assembly_bases}" >> $sample_log
    """
    }
}
process saveReferenceLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(reference_log)
    tuple val(sample_id),val(data_type),val(read_data),val(assembly_data),val(assembly_contigs),val(assembly_bases)

    script:

    if(!file(reference_log).isFile()){
        error "$reference_log doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${data_type}\t${read_data}\t${assembly_data}\t${assembly_contigs}\t${assembly_bases}" >> $reference_log
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
    echo "Sample_ID\tReference_ID\tPercent_Reference_Covered\tPercent_Query_Covered\tCategory\tYenta_SNPs\tgSNPs\tFiltered_Identity\tFiltered_Edge\tFiltered_Duplicated\tRejected_Density_1000\tRejected_Density_125\tRejected_Density_15" > ${results_directory}/MUmmer_DNADiff_Results.tsv
    echo -n "${results_directory}/MUmmer_DNADiff_Results.tsv"
    """
}
process saveDNADiffLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(diff_log)
    tuple val(sample_id),val(reference_id),val(percent_ref_cov),val(percent_query_cov),val(rsnps),val(category),val(gsnps),val(edge),val(iden),val(align),val(dup),val(het),val(density)

    script:

    if(!file(diff_log).isFile()){
        error "$diff_log doesn't exist..."
    } else{
    """
    echo "${sample_id}\t${reference_id}\t${percent_ref_cov}\t${percent_query_cov}\t${rsnps}\t${category}\t${gsnps}\t${edge}\t${iden}\t${align}\t${dup}\t${het}\t${density}" >> $diff_log
    """
    }
}
