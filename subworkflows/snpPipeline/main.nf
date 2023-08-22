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

// Check if parent folder exists
if(!output_directory.getParent().isDirectory()){
    error "Parent directory for output (--outroot) is not a valid directory [${output_directory.getParent()}]..."
} else if(!output_directory.isDirectory()){ // Check if output directory exists
    output_directory.mkdirs()
    mummer_directory.mkdirs()
    raw_mummer_directory.mkdirs()
} else{
    if(!mummer_directory.isDirectory()){ // Check if MUmmer directories exist
        mummer_directory.mkdirs()
        raw_mummer_directory.mkdirs()
    } else if(!raw_mummer_directory.isDirectory()){
        raw_mummer_directory.mkdirs()
    }
}

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

    snp_output = sample_pairwise
}

process mergeSNPs{

    output:
    stdout

    script:
    """
    echo "Isolate_ID\tData_Type\tRead_Data\tAssembly_Data\tAssembly_Contigs\tAssembly_Bases" > "${output_directory}/Query_Data.tsv"
    echo -n "${output_directory}/Query_Data.tsv"
    """
}
