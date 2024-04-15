// Subworkflow to run RefChooser for list of queries

// Set directory structure
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)

assembly_file = file("${log_directory}/Query_Assemblies.txt")
refchooser_matrix = file("${log_directory}/Refchooser_Matrix.tsv")
ref_id_file = file(params.ref_id_file)

workflow runRefChooser{
    take:
    query_data

    emit:
    reference_data

    main:

    // Run RefChooser
    refchooser_results = query_data
    .unique{it -> it[1]}
    .collect{it -> it[1]}
    | refChooser 
    | splitCsv | collect | flatten | collate(1) 
    | map{it -> tuple(it[0].toString(),null)}
   
    // Set reference_data channel 
    reference_data = query_data
    .map{it -> tuple(it[1].toString(),it[0])}
    .join(refchooser_results, by:0)
    .map{tuple(it[1],it[0])}
    .unique{it -> it[0]}.collect().flatten().collate(2)

    // Save reference data to file
    reference_data
    .collect{it -> it[0]}
    | saveRefIDs
}

process refChooser{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(assembly_paths)

    output:
    stdout

    script:

    ref_count = params.n_ref.toInteger()
    head_count = ref_count + 1
    assembly_file.write(assembly_paths.join('\n') + '\n')
    """
    $params.load_refchooser_module
    cd $log_directory
    refchooser metrics --sort Score $assembly_file sketch_dir > refchooser_results.txt
    refchooser maxtrix $assembly_file sketch_dir $refchooser_matrix

    column_data=\$(head -$head_count refchooser_results.txt | tail -$ref_count | cut -f7)
    head -$head_count refchooser_results.txt | tail -$ref_count | cut -f2,3,6,7,8 | grep -v Mean_Distance > refchooser.tsv
    if [[ \$(wc -l <<< "\$column_data") -gt 1 ]]; then
        echo "\$column_data" | paste -sd ',' -
    else
        echo -n "\$column_data"
    fi
    """
}

process saveRefIDs{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(ref_ids)

    script:
    ref_id_file.append(ref_ids.join('\n') + '\n')        
    """
    """
}