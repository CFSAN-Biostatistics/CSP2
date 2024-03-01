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

// Logging Processes//
// Takes a flattened list of snpdiffs files and generates a TSV report via python

process saveMUMmerLog{
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    input:
    val(mummer_data)

    output:
    stdout

    script:
 
    """
    echo "${mummer_data.join('\n')}"
    """
}