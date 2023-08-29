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
    reference_sample

    main:
    
    // Get reference isolate
    ref_path = sample_data | writeAssemblyPath | collect | flatten | first | refChooser
    
    sample_data.combine(ref_path).subscribe{println("Raw: $it")}
    reference_data = sample_data.combine(ref_path) | collect | flatten | collate(5) 
    | branch{
        
        same: "${it[3]}" == "/flash/storage/scratch/Robert.Literman/NextFlow/YENTA/GitHub/Yenta/NCBI_Clusters/Salmonella_PDS000043084.34/skesa_contigs/SRR10843746.fasta"
        return(tuple(it[0],it[1],it[2],it[3]))
        
        different: true
        return(it)}
    
    reference_data.same.subscribe{println("Same: $it")}
    reference_data.different.subscribe{println("Diff: $it")}
    reference_sample = reference_data.same
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