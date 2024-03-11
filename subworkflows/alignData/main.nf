// Subworkflow to run MUMmer for query/referece comparisons

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['align','assemble', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'align','assemble', 'screen', or 'snp', not ${params.runmode}..."
}

// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

// Set MUMmer directory structure
mummer_directory = file("${output_directory}/MUMmer_Output")
mum_coords_directory = file("${mummer_directory}/1coords")
mum_report_directory = file("${mummer_directory}/report")
mum_snps_directory = file("${mummer_directory}/snps")

// Set path to accessory scripts
mummerScript = file("$projectDir/bin/compileMUMmer.py")

// Set up modules if needed
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
params.load_mummer_module = params.mummer_module == "" ? "" : "module load -s ${params.mummer_module}"
params.load_bedtools_module = params.bedtools_module == "" ? "" : "module load -s ${params.bedtools_module}"

workflow alignGenomes{
    take:
    combined_data

    emit:
    return_mummer

    main:
    
    sample_pairwise = combined_data
    .filter{"${it[0]}" != "${it[2]}"} // Don't map things to themselves
    | runMUMmer | splitCsv
    
    return_mummer = sample_pairwise.collect().flatten().collate(3)
}

process runMUMmer{

    cpus = 1
    memory '4 GB'

    input:
    tuple val(query_name),val(query_fasta),val(ref_name),val(ref_fasta)

    output:
    stdout

    script:

    report_id = "${query_name}__vs__${ref_name}"

    // Ensure MUmmer directories exist
    if(!mummer_directory.isDirectory()){
        error "$mummer_directory does not exist..."
    } else{
        """
        module purge
        $params.load_mummer_module
        $params.load_python_module
        $params.load_bedtools_module

        cd ${mummer_directory}
        dnadiff -p ${report_id} ${ref_fasta} ${query_fasta}
        
        rm -rf ${mummer_directory}/${report_id}.mdelta
        rm -rf ${mummer_directory}/${report_id}.mcoords
        rm -rf ${mummer_directory}/${report_id}.1delta
        rm -rf ${mummer_directory}/${report_id}.delta
        rm -rf ${mummer_directory}/${report_id}.qdiff
        rm -rf ${mummer_directory}/${report_id}.rdiff
        rm -rf ${mummer_directory}/${report_id}.unref
        rm -rf ${mummer_directory}/${report_id}.unqry

        mv ${mummer_directory}/${report_id}.snps ${mum_snps_directory}
        mv ${mummer_directory}/${report_id}.report ${mum_report_directory}
        mv ${mummer_directory}/${report_id}.1coords ${mum_coords_directory}

        python ${mummerScript} "${query_name}" "${query_fasta}" "${ref_name}" "${ref_fasta}" "${output_directory}"       
        """
    }
}