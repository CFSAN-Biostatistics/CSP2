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

n_ref = params.n_ref.toInteger()
collate_num = n_ref+4

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
    hold_file = sample_data | writeAssemblyPath | collect | flatten | first 
    
    ref_path = refChooser(hold_file,n_ref) | splitCsv

    ref_ch = ref_path.map { line -> line.collect { it.trim() } } | view()
    reference_data = sample_data.filter { tuple -> ref_ch.any { refString -> tuple[3] == refString} } | view()
}

process refChooser{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(assembly_file)
    val(n_ref)

    output:
    stdout

    script:

    head_count = n_ref + 1

    """
    $params.load_refchooser_module
    cd $assembly_directory

    refchooser metrics --sort Score $assembly_file sketch_dir > refchooser_results.txt

    column_data=\$(head -$head_count refchooser_results.txt | tail -$n_ref | cut -f7)

    if [[ \$(wc -l <<< "\$column_data") -gt 1 ]]; then
        echo "\$column_data" | paste -sd ',' -
    else
        echo -n "\$column_data"
    fi
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