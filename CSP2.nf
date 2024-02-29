#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// CSP2 Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Assess run mode

// Runmode 'assemble'
//  - Requires --reads or --ref_reads
//  - Runs SKESA and summarzies output FASTA 

// Runmode 'align'
//  - Requires --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - Runs MUMmer, generates .snpdiffs, and alignment summary.
//      - If references are provided, all queries are aligned to all refs
//      - If no references are provided, query alignments are all-vs-all

// Runmode 'screen'
//  - Requires (1) --snpdiffs AND --ref_id or (2) --reads/--fasta AND --ref_reads/--ref_fasta
//  - Takes .snpdiffs files, applies QC, and generates SNP distance data between each query and reference pair

// Runmode 'snp'
//  - Requires (1) --snpdiffs AND --ref_id or (2) --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - If references are not provided, runs RefChooser to choose references (--n_ref sets how many references to choose)
//  - Takes .snpdiffs files, applies QC, and generates SNP distance data between all queries based on their alignment to each reference

if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Ensure data matches run mode
if ((run_mode == "assemble") && (params.reads == "" && params.ref_reads == "")) {
    error "Runmode is --assemble but no read data provided via --reads/--ref_reads"
} else if ((run_mode == "align") && (params.fasta == "" && params.reads == "")) {
    error "Runmode is --align but no query data provided via --fasta/--reads"
} else if ((run_mode == "screen") && (params.snpdiffs == "")){
    if((params.ref_fasta == "" && params.ref_reads == "")) {
        error "Runmode is --screen but no reference data provided via --snpdiffs/--ref_fasta/--ref_reads"
    } else if ((params.fasta == "" && params.reads == "")) {
        error "Runmode is --screen but no query data provided via --snpdiffs/--reads/--fasta"
    } 
} else if ((run_mode == "snp") && (params.snpdiffs == "")){
    if ((params.fasta == "" && params.reads == "")) {
        error "Runmode is --snp but no query data provided via --snpdiffs/--reads/--fasta"
    } 
} else{
    // Set directory structure
    if (params.outroot == "") {
        output_directory = file(params.out)
    } else {
        output_directory = file("${file(params.outroot)}/${params.out}")
    }
    
    if(output_directory.isDirectory()){
        error "${output_directory.getSimpleName()} (--out) already exists in ${output_directory.getParent()} (--outroot)..."
    } else if(!output_directory.getParent().isDirectory()){
        error "Parent directory for output (--outroot) is not a valid directory [${output_directory.getParent()}]..."
    } else{
        output_directory.mkdirs()
        if(run_mode != "assemble"){
            log_directory = file("${output_directory}/logs")
            log_directory.mkdirs()
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////






// Import modules
include {fetchData} from "./subworkflows/fetchData/main.nf"
include {saveIsolateLog} from "./subworkflows/logging/main.nf"
include {alignGenomes} from "./subworkflows/alignData/main.nf"
include {runScreen;runSNPPipeline} from "./subworkflows/snpdiffs/main.nf"
include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    // If run mode is 'assemble', tasks are complete
    if(run_mode != "assemble"){        
         if(run_mode == "align"){
            input_data.query_data.combine(input_data.reference_data) | alignGenomes // Align all queries against each reference and generate snpdiffs
        } else{
            // If run mode is 'screen' or 'snp' and references are provided, use them. If not, run RefChooser to generate references
            if(params.ref_reads == "" && params.ref_fasta == ""){ 
                if(params.snpdiffs == ""){
                    reference_data = runRefChooser(input_data.query_data) // Run RefChooser to generate references if none provided
                }else{
                    mummer_results = Channel.empty() // If snpdiffs are provided without other references, skip alignment
                }   
            } else{
                reference_data = input_data.reference_data
            }

            mummer_results = input_data.query_data.combine(reference_data) | alignGenomes // Align all queries against each reference and generate snpdiffs

            if(run_mode == "screen"){
                input_data.snpdiffs_data.concat(mummer_results) | runScreen // Compare snpdiffs to generate a summary
            }
            else if(run_mode == "snp"){
                input_data.snpdiffs_data.concat(mummer_results) | runSNPPipeline // Generate pairwise SNP distances and alignments against each reference
            }
        } 
    }
}