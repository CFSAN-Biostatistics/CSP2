#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Main script for running Yenta
// Params are read in from command line or from nextflow.config

// Ensure sample data is provided
params.reads = ""
params.fasta = ""
if(params.reads == "" && params.fasta == ""){
    error "Must provide sample isolate data (--reads and/or --fasta)"
}

// Ensure reference data is provided
params.ref_reads = ""
params.ref_fasta = ""
if(params.ref_reads == "" && params.ref_fasta == ""){
    error "Must provide reference isolate data (--ref_reads and/or --ref_fasta)"
}

// Create directory structure
params.outbase = "${projectDir}"
params.out = "YENTA_${new java.util.Date().getTime()}"
output_directory = file("${params.outbase}/${params.out}")
outroot = output_directory.getParent()

if(output_directory.isDirectory()){
    error "${output_directory} (--out) already exists..."
} else if(!outroot.isDirectory()){
    error "$outroot (--outbase) is not a valid directory..."
} else{
    output_directory.mkdirs()
    file("${output_directory}/Reference_Strain_Data").mkdirs()
    file("${output_directory}/Screening_Results").mkdirs()
}

// Import modules
include {fetchSampleData; fetchReferenceData} from "./subworkflows/fetchData/main.nf"
include {dnaDiff} from "./subworkflows/dnadiff/main.nf"

workflow{

    ////// 01: Collect paths to data and assemble read data if necessary ////// 
    sample_data = fetchSampleData()
    reference_data = fetchReferenceData() 

    ////// 02: Run MUmmer dnadiff on all sample x reference combos //////
    mummer_results = dnaDiff(sample_data,reference_data)
} 
