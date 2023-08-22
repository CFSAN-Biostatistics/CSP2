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
include {runSnpPipeline; runScreen } from "./subworkflows/dnadiff/main.nf"
include {compileSNPs} from "./subworkflows/snpPipeline/main.nf"

workflow{

    ////// 01: Collect paths to data and assemble read data if necessary ////// 
    sample_data = fetchSampleData()
    
    ////// 02: If reference data is provided, run the screening pipeline. Otherwise, run a SNP analysis //////
    if(params.ref_reads == "" && params.ref_fasta == ""){
        mummer_results = sample_data | collect | flatten | collate(4) | runSnpPipeline
        mummer_results.into { compileSNPs(it) }
    } else{
        reference_data = fetchReferenceData()
        runScreen(sample_data,reference_data)
    }
}
