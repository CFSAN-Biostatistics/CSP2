#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Main script for running Yenta
// Params are read in from command line or from nextflow.config

// Ensure sample data is provided
if(params.reads == "" && params.fasta == ""){
    error "Must provide isolate data for SNP analysis or screening (--reads and/or --fasta)"
}
// Set directory structure (--out for folder name; --outroot for parent dir)
if(params.outroot == ""){
    output_directory = file("${params.out}")
} else{
    output_directory = file("${file("${params.outroot}")}/${params.out}")
}

if(output_directory.isDirectory()){
    error "${output_directory.getSimpleName()} (--out) already exists in ${output_directory.getParent()} (--outroot)..."
} else if(!output_directory.getParent().isDirectory()){
    error "Parent directory for output (--outroot) is not a valid directory [${output_directory.getParent()}]..."
} else{
    output_directory.mkdirs()
    file("${output_directory}/Assemblies").mkdirs()
    file("${output_directory}/MUmmer_Output").mkdirs()
    file("${output_directory}/MUmmer_Output/Raw").mkdirs()
    file("${output_directory}/MUmmer_Output/Raw/Reports").mkdirs()
    file("${output_directory}/MUmmer_Output/Raw/1coords").mkdirs()
    file("${output_directory}/MUmmer_Output/Raw/SNPs").mkdirs()
    if(params.ref_reads == "" && params.ref_fasta == ""){
        file("${output_directory}/SNP_Analysis").mkdirs()
    }
}

// Import modules
include {fetchSampleData; fetchReferenceData} from "./subworkflows/fetchData/main.nf"
include {runSnpPipeline; runScreen; runAllvAll } from "./subworkflows/dnadiff/main.nf"
include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{

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
}
