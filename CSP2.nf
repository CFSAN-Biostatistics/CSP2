#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// CSP2 Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (!['align','assemble', 'screen', 'snp'].contains(params.runmode)){
    error "--runmode must be 'align','assemble', 'screen', or 'snp', not ${params.runmode}..."
}

// Ensure necessary data is provided given the run mode

// Runmode 'assemble'
//  - Requires: --reads/--ref_reads
//  - Runs SKESA and summarzies output FASTA 
if (params.runmode == "assemble"){
    if((params.reads == "") && (params.ref_reads == "")){
        error "Runmode is --assemble but no read data provided via --reads/--ref_reads"
    } 
}

// Runmode 'align'
//  - Requires: --reads/--fasta/--snpdiffs
//  - Optional: --ref_reads/--ref_fasta/--ref_id
//  - Runs MUMmer, generates .snpdiffs, and alignment summary.
//      - If references are provided via --ref_reads/--ref_fasta/--ref_id, non-reference samples are aligned to each reference
//      - If no references are provided, alignments are all-vs-all
//      - If --snpdiffs are provided, their FASTAs will be autodetected and, if present, used as queries or references as specified by --ref_reads/--ref_fasta/--ref_id
//      - Does NOT perform QC filtering

else if (params.runmode == "align"){
    if((params.fasta == "") && (params.reads == "") && (params.snpdiffs == "")){
        error "Runmode is --align but no query data provided via --fasta/--reads/--snpdiffs"
    } 
}

// Runmode 'screen'
//  - Requires: --reads/--fasta/--snpdiffs
//  - Optional: --ref_reads/--ref_fasta/--ref_id
//  - Generates .snpdiffs files (if needed), applies QC, and generates alignment summaries and SNP distance estimates
//      - If references are provided via --ref_reads/--ref_fasta/--ref_id, non-reference samples are aligned to each reference
//      - If no references are provided, alignments are all-vs-all
//      - If --snpdiffs are provided, (1) they will be QC filtered and included in the output report and (2) their FASTAs will be autodetected and, if present, used as queries or references as specified by --ref_reads/--ref_fasta/--ref_id

else if (params.runmode == "screen"){
    if((params.fasta == "") && (params.reads == "") && (params.snpdiffs == "")){
        error "Runmode is --screen but no query data provided via --snpdiffs/--reads/--fasta"
    }
}

// Runmode 'snp'
//  - Requires: --reads/--fasta/--snpdiffs
//  - Optional: --ref_reads/--ref_fasta/--ref_id
//  - If references are not provided, runs RefChooser using all FASTAs to choose references (--n_ref sets how many references to choose)
//  - Each query is aligned to each reference, and pairwise SNP distances for all queries are generated based on that reference
//  - Generates .snpdiffs files (if needed), applies QC, and generates SNP distance data between all queries based on their alignment to each reference
else if (params.runmode == "snp"){
    if((params.snpdiffs == "") && (params.fasta == "") && (params.reads == "")) {
        error "Runmode is --snp but no query data provided via --snpdiffs/--reads/--fasta"
    }
} 

// Set directory structure
if (params.outroot == "") {
    output_directory = file(params.out)
} else {
    out_root = file(params.outroot)
    output_directory = file("${out_root}/${params.out}")
}

// If the output directory exists, create a new subdirectory with the default output name ("CSP2_<TIME>")
if(!output_directory.getParent().isDirectory()){
    error "Parent directory for output (--outroot) is not a valid directory [${output_directory.getParent()}]..."
} else if(output_directory.isDirectory()){
    output_directory = file("${output_directory}/CSP2_${new java.util.Date().getTime()}")
    output_directory.mkdirs()
} else{
    output_directory.mkdirs()
}

if(params.tmp_dir != ""){
    tempdir = file(params.tmp_dir)

    if(tempdir.isDirectory()){
        temp_dir = file("${tempdir}/CSP2_${new java.util.Date().getTime()}_tmp")
        temp_dir.mkdirs()
        params.temp_dir = file(temp_dir)

    } else if(tempdir.getParent().isDirectory()){
        tempdir.mkdirs()
        params.temp_dir = tempdir
    } else{
        error "Parent directory for temp directory --tmp_dir (${params.tmp_dir}) is not a valid directory..."
    }
} else{
    params.temp_dir = ""
}

// Set MUMmer and SNP directories
mummer_directory = file("${output_directory}/MUMmer_Output")
snpdiffs_directory = file("${output_directory}/snpdiffs")
snp_directory = file("${output_directory}/SNP_Analysis")

// Set paths for output files 
isolate_data_file = file("${output_directory}/Isolate_Data.tsv")
screening_results_file = file("${output_directory}/Screening_Results.tsv")

// In --runmode assembly, results save to output_directory
if(params.runmode == "assemble"){
    ref_mode = false

    log_directory = file("${output_directory}")
    assembly_directory = file("${output_directory}")

    // Set dummy paths for log files
    screen_log_dir = file("${log_directory}/Screening_Logs")
    snp_log_dir = file("${log_directory}/SNP_Logs")
    ref_id_file = file("${log_directory}/Reference_IDs.txt")

} else{

    // Get reference mode
    if(params.ref_reads == "" && params.ref_fasta == "" && params.ref_id == ""){
        ref_mode = false
    } else{
        ref_mode = true
    }
    
    // Set directories
    log_directory = file("${output_directory}/logs")
    assembly_directory = file("${output_directory}/Assemblies")
    ref_id_file = file("${log_directory}/Reference_IDs.txt")

    // Create directories
    log_directory.mkdirs()
    mummer_directory.mkdirs()
    snpdiffs_directory.mkdirs()

    // Touch Reference_IDs.txt to establish it
    file(ref_id_file).text = ''

    // Set paths for log subdirectories
    screen_log_dir = file("${log_directory}/Screening_Logs")
    snp_log_dir = file("${log_directory}/SNP_Logs")
    
    // If --reads/--ref_reads are provided, prepare a directory for assemblies
    if((params.reads != "") || (params.ref_reads != "")){
        assembly_directory.mkdirs()
    }

    // If runmode is snp, prepare a directory for SNP analysis + logs
    if(params.runmode == "snp"){
        snp_directory.mkdirs()
        snp_log_dir.mkdirs()
    }

    // If runmode is screen, prepare a directory for screening logs
    if(params.runmode == "screen"){
        screen_log_dir.mkdirs()
    }
}

// Parameterize variables to pass between scripts
params.output_directory = file(output_directory)
params.log_directory = file(log_directory)
params.screen_log_dir = file(screen_log_dir)
params.snp_log_dir = file(snp_log_dir)
params.assembly_directory = file(assembly_directory)
params.mummer_directory = file(mummer_directory)
params.snpdiffs_directory = file(snpdiffs_directory)
params.snp_directory = file(snp_directory)
params.ref_id_file = file(ref_id_file)

params.ref_mode = ref_mode

// Set up modules if needed
params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
params.load_skesa_module = params.skesa_module == "" ? "" : "module load -s ${params.skesa_module}"
params.load_bedtools_module = params.bedtools_module == "" ? "" : "module load -s ${params.bedtools_module}"
params.load_bbtools_module = params.bbtools_module == "" ? "" : "module load -s ${params.bbtools_module}"
params.load_mummer_module = params.mummer_module == "" ? "" : "module load -s ${params.mummer_module}"
params.load_refchooser_module = params.refchooser_module == "" ? "" : "module load -s ${params.refchooser_module}"

// Save params to log file
params.each { key, value ->
    file("${log_directory}/CSP2_Params.txt") << "$key = $value\n"
}

//////////////////////////////////////////////////////////////////////////////////////////

// Import modules
include {fetchData} from "./subworkflows/fetchData/main.nf"
include {alignGenomes} from "./subworkflows/alignData/main.nf"
include {runScreen;runSNPPipeline} from "./subworkflows/snpdiffs/main.nf"
include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    query_data = input_data.query_data
    reference_data = input_data.reference_data
    snpdiffs_data = input_data.snpdiff_data

    // Create channel for pre-aligned data [(Query_1,Query_2),SNPDiffs_File]
    prealigned = snpdiffs_data
    .map { it -> tuple([it[0], it[1]].sort().join(',').toString(), it[2]) }
    .unique{it -> it[0]}
    
    // If run mode is 'assemble', tasks are complete
    if((params.runmode == "align") || (params.runmode == "screen") || (params.runmode == "snp")){

        // If there is no reference data, align all query_data against each other
        if(!ref_mode){

            if((params.runmode == "align") || (params.runmode == "screen")){
                seen_combinations = []
                
                to_align = query_data.combine(query_data) // Self-combine query data
                .collect().flatten().collate(4)
                .filter{it -> (it[1].toString() != "null") && (it[3].toString() != "null")} // Can't align without FASTA
                .filter{ it -> // Get unique combinations
        
                combination = ["${it[0]}", "${it[2]}"].sort()
                
                if(combination in seen_combinations) {
                    return false
                } else {
                    seen_combinations << combination
                    return true
                }}
            } 
            
            // If running SNP pipeline without references, run RefChooser to choose references
            else if(params.runmode == "snp"){
                reference_data = runRefChooser(query_data)
                
                to_align = query_data
                .combine(reference_data)
                .filter{it -> (it[1].toString() != "null") && (it[3].toString() != "null")} // Can't align without FASTA
            }
        }

        // If references are provided, align all queries against all references
        else{
            to_align = query_data
            .combine(reference_data)
            .filter{it -> (it[1].toString() != "null") && (it[3].toString() != "null")} // Can't align without FASTA
        }

        // Don't align things that are already aligned via --snpdiffs
        unaligned = to_align
        .map { it -> tuple([it[0], it[2]].sort().join(',').toString(),it[0], it[1], it[2], it[3]) }
        .unique{it -> it[0]}
        .join(prealigned, by:0, remainder:true)
        .filter{it -> it[5].toString() == "null"}
        .map{it -> [it[1], it[2], it[3], it[4]]}
        | collect | flatten | collate(4)

        all_snpdiffs = alignGenomes(unaligned,snpdiffs_data)
        .ifEmpty { error "No .snpdiffs to process..." }
        .collect().flatten().collate(3)
        
        if(params.runmode == "screen"){
            runScreen(all_snpdiffs)
        } else if(params.runmode == "snp"){
            runSNPPipeline(all_snpdiffs,reference_data)
        }
    }
}