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

    emit:
    reference_data

    main:
    
    // Get reference isolate
    ref_path = sample_data | writeAssemblyPath | collect | flatten | first | refChooser
    sample_data.map{it[3]}.subscribe{println("3: $it")}
    reference_data = sample_data.filter{"${it[3]}" == "${ref_path}"} | collect | flatten | collate(4)
}

process refChooser{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(assembly_file)

    output:
    stdout

    script:
    """
    $params.load_refchooser_module
    cd $assembly_directory
    refchooser metrics --sort Score $assembly_file sketch_dir > refchooser_results.txt && echo \$(head -2 refchooser_results.txt | tail -1 | cut -f7)
    """
}

process writeAssemblyPath{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    tuple val(sample_name),val(data_type),val(read_location),val(assembly_location)

    output:
    val(assembly_file)

    script:
    """
    echo "${assembly_location}" >> $assembly_file
    """
}