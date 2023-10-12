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
    hold_file = sample_data | map{it -> it[3]} | collect | writeAssemblyPath
    ref_path = refChooser(hold_file,n_ref) | splitCsv | collect | flatten | collate(1)
    
    combo_data = sample_data.combine(ref_path) | collect | flatten | collate(5)
    reference_data = combo_data.map{it -> tuple(it[0],it[1],it[2],it[3].toString(),it[4].toString())}.filter{it[3] == it[4]}.map{it->tuple(it[0],it[1],it[2],it[3])}
}

process writeAssemblyPath {
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(assemblies)

    output:
    stdout

    script:
    """
    echo "${assemblies.join('\n')}" > $assembly_file
    echo -n $assembly_file
    """
}

process refChooser{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(hold_file)
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