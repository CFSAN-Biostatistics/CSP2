// Logging Processes//

// Set file paths
snpdiffs_list_file = file(params.snpdiffs_list_file)
snpdiffs_summary_file = file(params.snpdiffs_summary_file)

// Set path to accessory scripts
saveSNPDiffs = file("$projectDir/bin/saveSNPDiffs.py")

process saveMUMmerLog{
// Takes: Flattened list of snpdiffs files and generates a TSV report via python
// Returns: Path to list of snpdiffs files
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(snpdiffs_paths)

    script:

    snpdiffs_list_file.write(snpdiffs_paths.join('\n'))
    """
    $params.load_python_module
    python ${saveSNPDiffs} "${snpdiffs_list_file}" "${snpdiffs_summary_file}"
    """
}