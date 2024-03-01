#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// CSP2 Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['assemble', 'align', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'assemble', 'align', 'screen', or 'snp', not ${params.runmode}..."
}

// Runmode 'assemble'
//  - Requires --reads or --ref_reads
//  - Runs SKESA and summarzies output FASTA 
if (run_mode == "assemble"){
    if((params.reads == "") && (params.ref_reads == "")){
        error "Runmode is --assemble but no read data provided via --reads/--ref_reads"
    } 
}

// Runmode 'align'
//  - Requires --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - Runs MUMmer, generates .snpdiffs, and alignment summary.
//      - If references are provided, all queries are aligned to all refs
//      - If no references are provided, query alignments are all-vs-all
else if (run_mode == "align"){
    if((params.fasta == "") && (params.reads == "")){
        error "Runmode is --align but no query data provided via --fasta/--reads"
    } 
}

// Runmode 'screen'
//  - Requires (1) --ref_id/--ref_reads/--ref_fasta AND (2) --reads/--fasta/--snpdiffs
//  - Takes .snpdiffs files, applies QC, and generates SNP distance data between each query and reference pair
else if (run_mode == "screen"){
    if((params.ref_fasta == "") && (params.ref_reads == "") && (params.ref_id == "")) {
        error "Runmode is --screen but no reference data provided via --ref_id/--ref_fasta/--ref_reads"
    } else if ((params.snpdiffs == "") && (params.fasta == "") && (params.reads == "")) {
        error "Runmode is --screen but no query data provided via --snpdiffs/--reads/--fasta"
    } else if(((params.snpdiffs != "") && (params.fasta == "") && (params.reads == "")) && (params.ref_id == "")) {
        error "Runmode is --screen, and only snpdiffs were provided, but no reference data provided via --ref_id"
    }
}

// Runmode 'snp'
//  - Requires (1) --snpdiffs AND --ref_id or (2) --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - If references are not provided, runs RefChooser to choose references (--n_ref sets how many references to choose)
//  - Takes .snpdiffs files, applies QC, and generates SNP distance data between all queries based on their alignment to each reference
else if (run_mode == "snp"){
    if((params.snpdiffs == "") && (params.fasta == "") && (params.reads == "")) {
        error "Runmode is --snp but no query data provided via --snpdiffs/--reads/--fasta"
    } else if(((params.snpdiffs != "") && (params.fasta == "") && (params.reads == "")) && (params.ref_id == "")) {
        error "Runmode is --snp, and only snpdiffs were provided, but no reference data provided via --ref_id"
    }
} 

// Runmode 'kmer' (IN DEVELOPMENT)
//  - Requires --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - Generates full kmer comparisons between all queries and references
//  - If no references are provided, query kmer comparisons are all-vs-all
    
else{
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
include {alignGenomes} from "./subworkflows/alignData/main.nf"
include {saveMUMmerLog} from "./subworkflows/logging/main.nf"

//include {saveIsolateLog} from "./subworkflows/logging/main.nf"
//include {runScreen;runSNPPipeline} from "./subworkflows/snpdiffs/main.nf"
//include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    // If run mode is 'assemble', tasks are complete
    if(run_mode != "assemble"){        
         
         if(run_mode == "align"){

            // If there is no reference data, align all query_data against each other
            if(params.ref_reads == "" && params.ref_fasta == "" && params.ref_id == ""){
                
                seen_combinations = []
                query_channel = input_data.query_data.combine(input_data.query_data).collect().flatten().collate(4)
                .filter{"${it[0]}" != "${it[2]}"} // Don't map things to themselves
                .filter { it ->
                combination = ["${it[0]}", "${it[2]}"].sort()
                if (combination in seen_combinations) {
                    return false
                    } else {
                        seen_combinations << combination
                        return true
                        }
                    }
                    mummer_results = query_channel | alignGenomes
            } else{
                mummer_results = input_data.query_data.combine(input_data.reference_data) | alignGenomes
            }
         } else if(run_mode == "screen"){
            mummer_results = input_data.query_data.combine(input_data.reference_data) | alignGenomes
            screen_results = input_data.snpdiffs_data.concat(mummer_results) | runScreen
        } else if(run_mode == "snp"){
            if(params.ref_reads == "" && params.ref_fasta == "" && params.ref_id == ""){
            
            }

        }
    }
}
        
        /*else{
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

*/