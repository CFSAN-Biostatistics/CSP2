// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}

assembly_directory = file("${output_directory}/Assemblies")

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

    emit:
    sample_data
    reference_data

    main:
    sample_data.map{it -> it[4]} | view | text | saveAs("${assembly_directory}/Assemblies.txt")
    ref_path = refChooser("${assembly_directory}/Assemblies.txt")

    reference_data = sample_data.branch{
        same: "${it[4]}" == "${ref_path}"
        return(it)}
}

// Log functions //
process refChooser{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    assembly_dir

    output:
    stdout

    script:
    """
    $params.load_refchooser_module
    cd $assembly_dir
    refchooser metrics --sort Score --threads 16 Assemblies.txt sketch_dir > refchooser_results.txt
    head -2 refchooser_results.txt | tail -1 | cut -f7 | echo
    """
}