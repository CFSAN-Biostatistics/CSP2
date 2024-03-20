// Subworkflow to run MUMmer for query/referece comparisons

// Set path variables
mummer_directory = file(params.mummer_directory)
snpdiffs_directory = file(params.snpdiffs_directory)

// Set path to accessory scripts
mummerScript = file("$projectDir/bin/compileMUMmer.py")

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
        $params.load_bbtools_module

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

        python ${mummerScript} "${query_name}" "${query_fasta}" "${ref_name}" "${ref_fasta}" "${mummer_directory}" "${snpdiffs_directory}"    
        """
    }
}