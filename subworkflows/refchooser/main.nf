// Subworkflow to run RefChooser for list of queries
// Params are passed from Yenta.nf or from the command line if run directly

// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

assembly_directory = file("${output_directory}/logs")
assembly_file = file("${assembly_directory}/Query_Assemblies.txt")

// Assess run mode
if (params.runmode == "") {
    if (params.snpdiffs != "") {
        run_mode = "screen" // If .snpdiffs are provided, generate a summary of query/reference alignments
    } else if((params.reads != "" || params.ref_reads != "") && (params.fasta == "" && params.ref_fasta == "")){
        run_mode = "assemble" // If only reads are provided, generate assemblies
    } else if(params.reads == "" && params.fasta == ""){
        error "No query data provided via --reads/--fasta/--snpdiffs" // Exit if no data is provided
    } else if (params.ref_fasta == "" && params.ref_reads == ""){
        run_mode = "snp" // If query data is provided without reference data, run the SNP pipeline with RefChooser
    } else if((params.reads != "" || params.fasta != "") && (params.ref_reads != "" || params.ref_fasta == "")){
        run_mode = "screen" // If query and reference data are provided, perform MUMmer alignment and generate a summary
    } else if((params.fasta == "" && params.reads == "") && (params.ref_fasta != "" || params.ref_reads != "")){
        error "Reference data provided via --ref_reads/--ref_fasta, but no query data provided by --reads/--fasta/--snpdiffs" // Exit if no query data is provided
    } 
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Get number of references
n_ref = "${params.n_ref}".toInteger()

// Set modules if necessary
params.load_refchooser_module = params.refchooser_module == "" ? "" : "module load -s ${params.refchooser_module}"

workflow runRefChooser{
    take:
    query_data

    emit:
    reference_data

    main:
    
    // Get reversed query
    reversed_query = query_data.map{it->tuple(it[1].toString(),it[0])}

    // Get reference isolate
    ref_path = refChooser(reversed_query.collect{it[0]},n_ref).splitCsv().collect().flatten().collate(1).map{tuple(it[0],null)}
    
    reference_data = ref_path
    .join(reversed_query, by:0)
    .map{tuple(it[2],it[0])}
    .collect().flatten().collate(2)
}

process refChooser{
    
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(assembly_paths)
    val(n_ref)

    output:
    stdout

    script:

    head_count = n_ref + 1

    """
    $params.load_refchooser_module
    cd $assembly_directory

    echo "${assembly_paths.join('\n')}" > $assembly_file
    refchooser metrics --sort Score $assembly_file sketch_dir > refchooser_results.txt

    column_data=\$(head -$head_count refchooser_results.txt | tail -$n_ref | cut -f7)
    head -$head_count refchooser_results.txt | tail -$n_ref | cut -f2,3,6,7,8 | grep -v Mean_Distance > refchooser.tsv
    if [[ \$(wc -l <<< "\$column_data") -gt 1 ]]; then
        echo "\$column_data" | paste -sd ',' -
    else
        echo -n "\$column_data"
    fi
    """
}