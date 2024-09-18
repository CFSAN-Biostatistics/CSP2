// Subworkflow to run RefChooser for list of queries

// Set directory structure
output_directory = file(params.output_directory)
log_directory = file(params.log_directory)
mash_directory = file(params.mash_directory)

workflow runRefChooser{
    take:
    query_data

    emit:
    reference_data

    main:

    // Make MASH sketches (1 CPU per query) and generate triangle (all CPUs)
    mash_refs = query_data
    .unique{it -> it[1]}
    .map { [ it[0], it[1] ] }
    | mashSkech 
    | collect
    | mashTriangle
    | chooseRefs
    | splitCsv | collect | flatten | collate(1) 

    reference_data = query_data
    .map{it -> tuple(it[1].toString(),it[0])}
    .join(mash_refs, by:0)
    .map{tuple(it[1],it[0])}
    .unique{it -> it[0]}.collect().flatten().collate(2)

    // Save reference data to file
    reference_data
    .collect{it -> it[0]}
    | saveRefIDs
}

process chooseRefs{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(mash_triangle)

    output:
    stdout

    script:

    ref_count = params.n_ref.toInteger()
    ref_script = file("${projectDir}/bin/chooseRefs.py")
    """
    $params.load_python_module  
    cd $mash_directory

    python $ref_script $ref_count $mash_triangle "${params.trim_name}"
    """
}

process mashTriangle{

    input:
    val(mash_sketches)

    output:
    stdout

    script:

    sketch_file = file("${mash_directory}/Mash_Sketches.txt")
    mash_triangle_file = file("${mash_directory}/Mash_Triangle")

    """
    $params.load_mash_module
    ls ${mash_directory}/*.msh > $sketch_file
    mash triangle -p ${params.cores} -l $sketch_file > $mash_triangle_file
    echo -n $mash_triangle_file
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
    echo -n "${mash_path}"
    """
}

process saveRefIDs{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(ref_ids)

    script:
    ref_id_file = file(params.ref_id_file)
    ref_id_file.append(ref_ids.join('\n') + '\n')        
    """
    """
}
