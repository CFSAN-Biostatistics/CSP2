// Screening and SNP Pipeline processing
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)
screen_log_dir = file(params.screen_log_dir)
snp_log_dir = file(params.snp_log_dir)
snp_directory = file(params.snp_directory)

if(params.tmp_dir == ""){
    temp_dir = ""
} else{
    temp_dir = file(params.temp_dir)
}
ref_id_file = file(params.ref_id_file)

ref_mode = params.ref_mode

// Assess whether to rescue edge-filtered SNPs
edge_rescue = "${params.rescue}" == "norescue" ? "norescue" : "rescue"

// Set paths for output files 
all_snpdiffs_list = file("${log_directory}/All_SNPDiffs.txt")
snp_dirs_list = file("${log_directory}/SNP_Dirs.txt")
screening_results_file = file("${output_directory}/Screening_Results.tsv")
isolate_data_file = file("${output_directory}/Isolate_Data.tsv")
snpdiffs_summary_file = file("${output_directory}/Raw_MUMmer_Summary.tsv")

// Get QC thresholds
min_cov = params.min_cov.toFloat()
min_length = params.min_len.toInteger()
min_iden = params.min_iden.toFloat()
reference_edge = params.ref_edge.toInteger()
query_edge = params.query_edge.toInteger()
max_missing = params.max_missing.toFloat()
n_ref = params.n_ref.toInteger()

// iqTree bootstraps
np_boot = params.b.toInteger()
uf_boot = params.bb.toInteger()
boot_sum = np_boot + uf_boot

workflow {
    main:
    // Run SNP pipeline
    runSNPPipeline(query_data: all_snpdiffs, reference_data: ref_id_file)
}

workflow runScreen {
    
    take:
    all_snpdiffs

    main:

    all_snpdiffs
    .unique{it -> it[2]} 
    .collect() 
    | screenSNPDiffs
}

process screenSNPDiffs{

    input:
    val(all_snpdiffs)

    script:

    screenDiffs = file("${projectDir}/bin/screenSNPDiffs.py")
    """
    $params.load_python_module
    $params.load_bedtools_module
    python $screenDiffs --snpdiffs_file "${all_snpdiffs_list}" --log_dir "${screen_log_dir}" --min_cov "${min_cov}" --min_len "${min_length}" --min_iden "${min_iden}" --ref_edge "${reference_edge}" --query_edge "${query_edge}" --density_windows "${params.dwin}" --max_snps "${params.wsnps}" --trim_name "${params.trim_name}" --output_file "${screening_results_file}" --ref_id "${ref_id_file}" --tmp_dir "${temp_dir}"
    """
} 

workflow runSNPPipeline{
    take:
    all_snpdiffs
    reference_data

    main:
    
    query_snpdiffs = all_snpdiffs.map{tuple(it[0],it[2])}
    ref_snpdiffs = all_snpdiffs.map{tuple(it[1],it[2])}

    stacked_snpdiffs = query_snpdiffs.concat(ref_snpdiffs)
    .collect().flatten().collate(2)

    snp_dirs = stacked_snpdiffs
    .combine(reference_data)
    .filter{it -> it[0].toString() == it[2].toString()}
    .map{it -> tuple(it[0],it[1])}
    .groupTuple(by:0)
    .map { ref, diff_files -> tuple( ref.toString(), diff_files.collect() ) }
    | runSnpPipeline

    snp_trees = "${params.notree}" != "notree" ? Channel.empty() : snp_dirs | runiqtree
}

process compileResults{

    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(snp_directories)

    script:

    compile_script = file("${projectDir}/bin/compileSNPResults.py")
    snp_dirs_list.write(snp_directories.join("\n")+ "\n")
    """
    $params.load_python_module
    python $compile_script --snp_dirs_file "${snp_dirs_list}" --output_directory "${snp_directory}" --isolate_data_file "${isolate_data_file}" --mummer_data_file "${snpdiffs_summary_file}"
    """
}

process runSnpPipeline{

    input:
    tuple val(reference_id),val(diff_files)

    output:
    stdout

    script:

    snp_script = file("${projectDir}/bin/runSNPPipeline.py")

    // Set + create output directory
    snp_dir = file("${snp_directory}/${reference_id}")
    snp_dir.mkdirs()

    // Write SNPDiffs list
    out_snpdiffs = file("${snp_dir}/SNPDiffs.txt")
    out_snpdiffs.write(diff_files.join("\n")+ "\n")
    """
    $params.load_python_module
    $params.load_bedtools_module
    python $snp_script --reference_id "${reference_id}" --output_directory "${snp_dir}" --snpdiffs_file "${out_snpdiffs}" --log_directory "${snp_log_dir}" --min_cov "${min_cov}" --min_len "${min_length}" --min_iden "${min_iden}" --ref_edge "${reference_edge}" --query_edge "${query_edge}" --density_windows "${params.dwin}" --max_snps "${params.wsnps}" --trim_name "${params.trim_name}" --max_missing "${max_missing}" --tmp_dir "${temp_dir}" --rescue "${edge_rescue}"
    echo -n $snp_dir
    """
}

process runiqtree{

    input:
    val(snp_dir)

    output:
    stdout

    script:
    preserved_alignment = file("${snp_dir}/snpma_preserved.fasta")
    raw_tree_file = file("${snp_dir}/snpma_preserved.fasta.treefile")
    tree_file = boot_sum == 0 ? file("${snp_dir}/snpma_preserved.noboot.tre") : file("${snp_dir}/snpma_preserved.best.withboot.tre")

    // Ensure alignment exist
    if(!preserved_alignment.isFile()){
        error "$preserved_alignment does not exist..."
    } else{
        if(boot_sum == 0){
            bootstrap = ""
        } else{
            bootstrap = np_boot != 0 ? "-b ${np_boot}":"-bb ${uf_boot}"
        }
        """
        iqtree -T AUTO -m $params.model $bootstrap -s $preserved_alignment
        ln -s $raw_tree_file $tree_file
        echo -n $tree_file
        """
    }
}   