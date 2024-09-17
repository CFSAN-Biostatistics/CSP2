// Subworkflow to run RefChooser for list of queries

// Set directory structure
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)
mash_directory = file(params.mash_directory)

assembly_file = file("${log_directory}/Query_Assemblies.txt")
sketch_file = file("${mash_directory}/Mash_Sketches.txt")
mash_triangle = file("${mash_directory}/Mash_Triangle")

ref_id_file = file(params.ref_id_file)

ref_script = file("$projectDir/bin/chooseRefs.py")

workflow runNewRefChooser{
    take:
    query_data

    emit:
    reference_data

    main:

    // Run MASH
    mash_triangle = query_data
    .unique{it -> it[1]}
    | mashSkech    
    .collect{it -> it[0]}
    | saveMashPaths
    | mashTriangle
}

process mashTriangle{

    input:
    val(sketch_file)

    output:
    stdout

    script:
    """
    $params.load_mash_module
    mash triangle -p ${params.cores} -l $sketch_file > $mash_triangle
    echo $mash_triangle
    """
}

process mashSkech{
    cpus = 1

    input:
    tuple val(query_name),val(query_fasta)
    
    output:
    stdout

    script:

    mash_path = "${mash_directory}/${query_name}.msh"
    """
    $params.load_mash_module
    mash sketch -s 10000 -p 1 -o $mash_path $query_fasta
    echo "${mash_path}"
    """
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
    assembly_file.write(assembly_paths.join('\n') + '\n')
    """
    $params.load_refchooser_module
    cd $log_directory
    
    refchooser metrics --sort Score $assembly_file ${log_directory}/sketch_dir > ${log_directory}/refchooser_results.txt
    refchooser matrix $assembly_file ${log_directory}/sketch_dir ${log_directory}/refchooser_matrix.txt
    
    $params.unload_refchooser_module
    $params.load_python_module
    python $ref_script $ref_count refchooser_results.txt refchooser_matrix.txt
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

process saveMashPaths{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(mash_paths)

    script:
    sketch_file.append(mash_paths.join('\n') + '\n')        
    """
    echo "${sketch_file}"
    """
}