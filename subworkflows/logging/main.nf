// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Set directory structure
if(params.outroot == "") {
    output_directory = file(params.out)
} else {
    output_directory = file("${file(params.outroot)}/${params.out}")
}

log_directory = file("${output_directory}/logs")
snpdiffs_directory = file("${output_directory}/snpdiffs")

mummer_directory = file("${output_directory}/MUMmer_Output")
mum_coords_directory = file("${mummer_directory}/1coords")
mum_report_directory = file("${mummer_directory}/report")
mum_snps_directory = file("${mummer_directory}/snps")

snpdiffs_list_file = file("${log_directory}/All_SNPDiffs.txt")
snpdiffs_summary_file = file("${output_directory}/Raw_Alignment_Summary.tsv")

// Set path to accessory scripts
saveSNPDiffs = file("$projectDir/bin/saveSNPDiffs.py")

// Load modules
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"

// Logging Processes//
// Takes a flattened list of snpdiffs files and generates a TSV report via python

process saveMUMmerLog{
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(snpdiffs_paths)

    output:

    script:

    """
    $params.load_python_module
    cd $mummer_directory

    echo "${snpdiffs_paths.join('\n')}" > $snpdiffs_list_file
    python ${saveSNPDiffs} "${snpdiffs_list_file}" "${snpdiffs_summary_file}"
    """
}