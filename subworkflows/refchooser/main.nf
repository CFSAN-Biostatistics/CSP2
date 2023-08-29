// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}

assembly_directory = file("${output_directory}/Assemblies")
assembly_file = file("${output_directory}/Assemblies/Assemblies.txt")

// Set modules if necessary
params.refchooser_module = ""
if(params.refchooser_module == ""){
    params.load_refchooser_module = ""
} else{
    params.load_refchooser_module = "module load -s ${params.refchooser_module}"
}

workflow runRefChooser{
    take:
    sample_data

    //emit:
    //sample_data
    //reference_data

    main:
    
    // Create assembly list
    ref_path = refChooser(assembly_file)

    reference_data = sample_data.branch{
        same: "${it[4]}" == "${ref_path}"
        return(it)}
}

process refChooser{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    assembly_file

    output:
    stdout

    script:
    """
    $params.load_refchooser_module
    cd $assembly_directory
    refchooser metrics --sort Score $assembly_file sketch_dir > refchooser_results.txt
    head -2 refchooser_results.txt | tail -1 | cut -f7 | echo
    """
}