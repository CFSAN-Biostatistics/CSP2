// Subworkflow to run MUmmer for query/referece comparisons
// Params are passed from Yenta.nf or from the command line if run directly

// Set output paths
params.out = "Yenta_${new java.util.Date().getTime()}"
params.outroot = ""
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}

mummer_directory = file("${output_directory}/MUmmer_Output")
snp_directory = file("${output_directory}/SNP_Analysis")

// Set path to snp script
snp_script = file("$projectDir/bin/mergeSNPs.py")

params.python_module = ""
if(params.python_module == ""){
    params.load_python_module = ""
} else{
    params.load_python_module = "module load -s ${params.python_module}"
}

workflow compileSNPs{
    take:
    sample_pairwise

    emit:
    snp_output

    main:

    snp_output = mergeSNPs()
}

process mergeSNPs{

    output:
    snp_directory

    script:
    """
    ${params.load_python_module}
    python $snp_script $mummer_directory $snp_directory
    """
}
