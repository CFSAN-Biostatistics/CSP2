#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// CSP2 Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (['align','assemble', 'filter', 'screen', 'snp'].contains(params.runmode)) {
    run_mode = "${params.runmode}"
} else {
    error "--runmode must be 'align','assemble', 'filter', 'screen', or 'snp', not ${params.runmode}..."
}

// Ensure necessary data is provided given the run mode

// Runmode 'assemble'
//  - Requires: --reads/--ref_reads
//  - Runs SKESA and summarzies output FASTA 
if (run_mode == "assemble"){
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

else if (run_mode == "align"){
    if((params.fasta == "") && (params.reads == "") && (params.snpdiffs == "")){
        error "Runmode is --align but no query data provided via --fasta/--reads/--snpdiffs"
    } 
}

// Runmode 'filter'
//  - Requires: --snpdiffs
// Runs QC on --snpdiffs to generate SNP distance data
else if (run_mode == "filter"){
    if(params.snpdiffs == ""){
        error "Runmode is --filter but no snpdiffs provided via --snpdiffs"
    } 
}

// Runmode 'screen'
//  - Requires: --reads/--fasta/--snpdiffs
//  - Optional: --ref_reads/--ref_fasta/--ref_id
//  - Generates .snpdiffs files (if needed), applies QC, and generates alignment summaries and SNP distance estimates
//      - If references are provided via --ref_reads/--ref_fasta/--ref_id, non-reference samples are aligned to each reference
//      - If no references are provided, alignments are all-vs-all
//      - If --snpdiffs are provided, (1) they will be QC filtered and included in the output report and (2) their FASTAs will be autodetected and, if present, used as queries or references as specified by --ref_reads/--ref_fasta/--ref_id

else if (run_mode == "screen"){
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
else if (run_mode == "snp"){
    if((params.snpdiffs == "") && (params.fasta == "") && (params.reads == "")) {
        error "Runmode is --snp but no query data provided via --snpdiffs/--reads/--fasta"
    }
} 

// Runmode 'kmer' (IN DEVELOPMENT)
//  - Requires --reads/--fasta
//  - Optional: --ref_reads/--ref_fasta
//  - Generates full kmer comparisons between all queries and references
//  - If no references are provided, query kmer comparisons are all-vs-all
    
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
    
    // In --runmode assembly, results save to output_directory
    if(run_mode != "assemble"){
        log_directory = file("${output_directory}/logs")
        log_directory.mkdirs()
        
        // If --reads/--ref_reads are provided, make a directory for assemblies
        if((params.reads != "") || (params.ref_reads != "")){
            assembly_directory = file("${output_directory}/Assemblies")
            assembly_directory.mkdirs()    
        }

        // If --fasta/--reads/--snpdiffs are provided, make directories for alignment data
        if((params.fasta != "") || (params.reads != "") || (params.snpdiffs != "")){
            mummer_directory = file("${output_directory}/MUMmer_Output")
            mum_coords_directory = file("${mummer_directory}/1coords")
            mum_report_directory = file("${mummer_directory}/report")
            mum_snps_directory = file("${mummer_directory}/snps")   
            snpdiffs_directory = file("${output_directory}/snpdiffs") 
            mummer_directory.mkdirs()
            mum_coords_directory.mkdirs()
            mum_report_directory.mkdirs()
            mum_snps_directory.mkdirs()
            snpdiffs_directory.mkdirs()
        }
        
        if(run_mode == "snp"){
            snp_directory = file("${output_directory}/SNP_Analysis")
            snp_directory.mkdirs()
        }
    }  
}

//////////////////////////////////////////////////////////////////////////////////////////

// Import modules
include {fetchData} from "./subworkflows/fetchData/main.nf"
include {alignGenomes} from "./subworkflows/alignData/main.nf"
include {saveMUMmerLog} from "./subworkflows/logging/main.nf"
include {runScreen} from "./subworkflows/snpdiffs/main.nf"

//include {saveIsolateLog} from "./subworkflows/logging/main.nf"
//include {runScreen;runSNPPipeline} from "./subworkflows/snpdiffs/main.nf"
//include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    // Get pre-aligned data from snpdiffs
    input_data.snpdiff_data.map { it -> ["${it[0]}", "${it[2]}"].sort() }.unique().toList().set { sortedUniqueList }  
    
    // If run mode is 'assemble', tasks are complete
    if(run_mode == "align"){

        // If there is no reference data, align all query_data against each other
        if(params.ref_reads == "" && params.ref_fasta == "" && params.ref_id == ""){
            
            seen_combinations = []
            
            to_align = input_data.query_data
            .combine(input_data.query_data).collect().flatten().collate(4)
            .filter { it ->
    
            combination = ["${it[0]}", "${it[2]}"].sort()
            
            if (combination in seen_combinations) {
                return false
            } else {
                seen_combinations << combination
                return true
            }}
        } else{
            to_align = input_data.query_data.combine(input_data.reference_data)
        }
        
        mummer_results = to_align | alignGenomes | collect | flatten | collate(3)
        saveMUMmerLog(mummer_results.collect{it[2]})

    } else if(run_mode == "screen"){

          
            
            snpdiff_aligned = 
            .map { it ->
            combination = ["${it[0]}", "${it[2]}"].sort()
            
            if (combination in pre_aligned) {
                return false
            } else {
                pre_aligned << combination
                return false
            }
            }
        } else{
            to_align = input_data.query_data.combine(input_data.reference_data)
        }
        mummer_results = input_data.query_data.combine(input_data.reference_data) | alignGenomes | collect | flatten | collate(3)
        all_snpdiffs = input_data.snpdiff_data.concat(mummer_results).collect().flatten().collate(3)
        saveMUMmerLog(all_snpdiffs.collect{it[2]})

        runScreen(all_snpdiffs,input_data.reference_data)
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