#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Main script for running Yenta
// Params are read in from command line or from nextflow.config

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


// Ensure data before continuing
if ((run_mode == "assemble") && (params.reads == "" && params.ref_reads == "")) {
    error "Runmode is --assemble but no read data provided via --reads/--ref_reads"
} else if ((run_mode == "align") && (params.fasta == "" && params.reads == "")) {
    error "Runmode is --align but no query data provided via --fasta/--reads"
}  else if ((run_mode == "align") && (params.ref_fasta == "" && params.ref_reads == "")) {
    error "Runmode is --align but no reference data provided via --ref_fasta/--ref_reads"
} else if ((run_mode == "screen") && (params.snpdiffs == "" && params.fasta == "" && params.reads == "")) {
    error "Runmode is --screen but no query data provided via --fasta/--reads/--snpdiffs"
}  else if ((run_mode == "screen") && (params.snpdiffs == "" && params.ref_fasta == "" && params.ref_reads == "")) {
    error "Runmode is --screen but no reference data provided via --ref_fasta/--ref_reads/--snpdiffs"
} else if ((run_mode == "snp") && (params.snpdiffs == "" && params.fasta == "" && params.reads == "")) {
    error "Runmode is --snp but no query data provided via --snpdiffs/--fasta/--reads"
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

// Import modules
include {fetchData} from "./subworkflows/fetchData/main.nf"
include {saveIsolateLog} from "./subworkflows/logging/main.nf"
include {alignGenomes} from "./subworkflows/alignData/main.nf"
//include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    // If run mode is 'assemble', tasks are complete
    if(run_mode != "assemble"){
        
        // Get query and reference data if provided
        query_data = input_data.query_data
        reference_data = input_data.reference_data

        if(run_mode == "align"){
            mummer_results = query_data.combine(reference_data) | alignGenomes // Align all queries against each reference and generate snpdiffs
            mummer_results.map{it -> it.join("\t")}.collect() | saveIsolateLog // Save isolate data
        } 
        
        else{
            
            // Fetch snpdiffs input if provided
            all_snp_diffs = input_data.snpdiffs_data.concat(mummer_results)

            if(run_mode == "screen"){
                mummer_results = query_data.combine(reference_data) | alignGenomes // Align all queries against each reference and generate snpdiffs
            } 
            
            else if(run_mode == "snp"){
                mummer_results = query_data.combine(reference_data) | alignGenomes // Align all queries against each reference and generate snpdiffs
            }
        
        
        } 
    }

}




"""



}

    ////// Read in sample data ///////
    sample_data = fetchSampleData()
     
    ////// Check run mode //////
    if(params.ref_reads != "" || params.ref_fasta != ""){   // If --ref_reads/--ref_fasta are set, run in reference screener mode
        reference_data = fetchReferenceData(params.ref_reads,params.ref_fasta)
        runScreen(sample_data,reference_data)
    } else{ // SNP Pipeline mode
        
        if(params.snp_ref_reads != "" || params.snp_ref_fasta != ""){    // If --snp_ref_reads/--snp_ref_fasta are set, fetch user-selected references
            reference_data = fetchReferenceData(params.snp_ref_reads,params.snp_ref_fasta)} 
        else{ // Get reference from RefChooser
            reference_data = runRefChooser(sample_data)}

        runSnpPipeline(sample_data,reference_data)
    }
"""