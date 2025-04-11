#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// CSP2 Main Script
// Params are read in from command line or from nextflow.config and/or conf/profiles.config

// Check if help flag was passed
help1 = "${params.help}" == "nohelp" ? "nohelp" : "help"
help2 = "${params.h}" == "nohelp" ? "nohelp" : "help"

def printHelp() {
    println """
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CSP2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Global default params:

  --out            Set name for output folder/file prefixes (Default: CSP2_<timestamp>)
  --outroot        Set output parent directory (Default: CWD; Useful to hardset in nextflow.config if 
                   you want all output go to the same parent folder, with unique IDs set by --out)
  --tmp_dir        Manually specify a TMP directory for pybedtools output
  --help/--h       Display this help menu


CSP2 can run in the following run modes:

  --runmode        Run mode for CSP2:

                   - assemble: Assemble read data (--reads/--ref_reads) into FASTA using SKESA
  
                   - align: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta), 
                            run MUMmer alignment analysis for each query/ref combination
  
                   - screen: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta) 
                             and/or MUMmer output (.snpdiffs), create a report for raw SNP 
                             distances between each query and reference assembly
  
                   - snp: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta) 
                          and/or MUMmer output (.snpdiffs), generate alignments and pairwise 
                          distances for all queries based on each reference dataset

                   - locus: Given query data (--reads/--fasta) and reference loci (--ref_fasta) 
                          and/or MUMmer output (.snpdiffs), detect loci in each query and, if present,
                          extract and align them
Input Data:

  --fasta          Location for query isolate assembly data (.fasta/.fa/.fna). Can be a list of files, a path 
                   to a signle single FASTA, or a path to a directories with assemblies. 
  --ref_fasta      Location for reference isolate assembly data (.fasta/.fa/.fna). Can be a list of files, a 
                   path to a signle single FASTA, or a path to a directories with assemblies. 
  --reads          Directory or list of directories containing query isolate read data
  --readext        Read file extension (Default: fastq.gz)
  --forward        Forward read file suffix (Default: _1.fastq.gz)
  --reverse        Reverse read file suffix (Default: _2.fastq.gz)
  
  --ref_reads      Directory or list of directories containing reference isolate read data
  --ref_readext    Reference read file extension (Default: fastq.gz)
  --ref_forward    Reference forward read file suffix (Default: _1.fastq.gz)
  --ref_reverse    Reference reverse read file suffix (Default: _2.fastq.gz)

  --snpdiffs       Location for pre-generated snpdiffs files (List of snpdiffs files, directory with snpdiffs)

  --ref_id         IDs to specify reference sequences (Comma-separated list; e.g., Sample_A,Sample_B,Sample_C)

  --trim_name      A common string to remove from all sample IDs (Default: ''; Useful if all assemblies end in 
                   something like "_contigs_skesa.fasta")

  --n_ref          If running in --runmode snp, the number of reference genomes for CSP2 to select if none are provided (Default: 1)

  --exclude        A comma-separated list of IDs to remove prior to analysis (Useful for removing low quality 
                   isolates in combination with --snpdiffs)

QC variables:

  --min_cov        Only consider queries if the reference genome is covered by at least <min_cov>% (Default: 85)
  --min_len        Only consider SNPs from contig alignments longer than <min_len> bp (Default: 500)
  --min_iden       Only consider SNPs from alignments with at least <min_iden> percent identity (Default: 99)
  --dwin           A comma-separated set of window sizes for SNP density filters (Default: 1000,125,15; Set --dwin 0 to disable density filtering)
  --wsnps          A comma-separated set of maximum SNP counts per window above (Default: 3,2,1)
  --max_missing    If running in --runmode snp, mask SNPs where data is missing or purged from <max_missing>% of isolates (Default: 50)

Edge Trimming:

  --ref_edge       Don't include SNPs that fall within <ref_edge>bp of a reference contig edge (Default: 150) 
  --query_edge     Don't include SNPs that fall within <query_edge>bp of a query contig edge (Default: 150)
  --rescue         If flagged (Default: not flagged), sites that were filtered out due solely to query edge proximity are rescued if 
                   the same reference position is covered more centrally by another query
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example Commands: 

1) Run CSP2 in SNP Pipeline mode using all the FASTA from /my/data/dir, and choose 3 references

nextflow run CSP2.nf --runmode snp --fasta /my/data/dir --n_ref 3

2) Screen all the paired-end .fastq files from /my/read/dir against the reference isolate in /my/reference/isolates.txt

nextflow run CSP2.nf --runmode screen --ref_fasta /my/reference/isolates.txt --reads /my/read/dir --readext .fastq --forward _1.fastq --reverse _2.fastq

3) Re-run the SNP pipeline using old snpdiffs files after changing the density filters and removing a bad sample

nextflow run CSP2.nf --runmode snp --snpdiffs /my/old/analysis/snpdiffs --dwin 5000,2500,1000 --wsnps 6,4,2 --ref_id Sample_A --exclude Sample_Q --out HQ_Density 

4) Run in assembly mode and use HPC modules specified in profiles.config (NOTE: Setting the profile in nextflow uses a single hyphen (-) as compared to other arguments (--))

nextflow run CSP2.nf -profile myHPC --runmode assemble --reads /my/read/dir --out Assemblies

5) Run in SNP pipeline mode using SLURM and use the built in conda environment (NOTE: For local jobs using conda, use -profile standard_conda)

nextflow run CSP2.nf -profile slurm_conda --runmode snp --fasta /my/data/dir --out CSP2_Conda

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
    System.exit(0)
}

if (help1 == "help") {
    printHelp()
} else if(help2 =="help"){
    printHelp()
}

// Assess run mode
if (params.runmode == "") {
    error "--runmode must be specified..."
} else if (!['align','assemble', 'screen', 'snp','locus','conda_init'].contains(params.runmode)){
    error "--runmode must be 'align','assemble', 'screen','locus', or 'snp', not ${params.runmode}..."
}

// If runmode is conda_init, launch a local process to spurn the generation of the conda environment and exit
if (params.runmode != "conda_init") {

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

    // Runmode 'locus'
    //  - Requires: --reads/--fasta/--snpdiffs/--ref_fasta
    //  - Generates .snpdiffs files (if needed), applies QC, and aligns loci extracted from queries if present
    else if (params.runmode == "locus"){
        if(params.snpdiffs == ""){
            if((params.fasta == "") && (params.reads == "")) {
                error "Runmode is --locus but no query data provided via --snpdiffs/--reads/--fasta"
            } else if(params.ref_fasta == ""){
                error "Runmode is --locus but no locus data provided via --snpdiffs/--ref_fasta"
            }
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
        mummer_log_directory = file("${log_directory}/MUMmer_Logs")
        mash_dir = file("${log_directory}/sketch_dir")

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
        mummer_log_directory = file("${log_directory}/MUMmer_Logs")
        ref_id_file = file("${log_directory}/Reference_IDs.txt")

        // Create directories
        log_directory.mkdirs()
        mummer_directory.mkdirs()
        mummer_log_directory.mkdirs()
        snpdiffs_directory.mkdirs()

        // Touch Reference_IDs.txt to establish it
        file(ref_id_file).text = ''

        // Set paths for log subdirectories
        screen_log_dir = file("${log_directory}/Screening_Logs")
        snp_log_dir = file("${log_directory}/SNP_Logs")
        mash_dir = file("${log_directory}/sketch_dir")

        // If --reads/--ref_reads are provided, prepare a directory for assemblies
        if((params.reads != "") || (params.ref_reads != "")){
            assembly_directory.mkdirs()
        }

        // If runmode is snp, prepare a directory for SNP analysis + logs
        if(params.runmode == "snp"){
            snp_directory.mkdirs()
            snp_log_dir.mkdirs()
            if(!ref_mode){
                mash_dir.mkdirs()
            }        
        }

        // If runmode is locus, prepare a directory for SNP analysis + logs
        if(params.runmode == "locus"){
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
    params.mummer_log_directory = file(mummer_log_directory)
    params.snpdiffs_directory = file(snpdiffs_directory)
    params.snp_directory = file(snp_directory)
    params.ref_id_file = file(ref_id_file)
    params.mash_directory = file(mash_dir)

    params.ref_mode = ref_mode

    // Set up modules if needed
    params.load_python_module = params.python_module == "" ? "" : "module load -s ${params.python_module}"
    params.load_skesa_module = params.skesa_module == "" ? "" : "module load -s ${params.skesa_module}"
    params.load_bedtools_module = params.bedtools_module == "" ? "" : "module load -s ${params.bedtools_module}"
    params.load_bbtools_module = params.bbtools_module == "" ? "" : "module load -s ${params.bbtools_module}"
    params.load_mummer_module = params.mummer_module == "" ? "" : "module load -s ${params.mummer_module}"
    params.load_mash_module = params.mash_module == "" ? "" : "module load -s ${params.mash_module}"

    // Save params to log file
    params.each { key, value ->
        file("${log_directory}/CSP2_Params.txt") << "$key = $value\n"
    }
} else{
    params.output_directory = "./"
    params.log_directory = "./"
    params.screen_log_dir = "./"
    params.snp_log_dir = "./"
    params.assembly_directory = "./"
    params.mummer_directory = "./"
    params.mummer_log_directory = "./"
    params.snpdiffs_directory = "./"
    params.snp_directory = "./"
    params.ref_id_file = "./"
    params.mash_directory = "./"
    params.ref_mode = false
    params.load_python_module = ""
    params.load_skesa_module = ""
    params.load_bedtools_module = ""
    params.load_bbtools_module = ""
    params.load_mummer_module = ""
    params.load_mash_module = ""
}

//////////////////////////////////////////////////////////////////////////////////////////

// Import modules
include {fetchData} from "./subworkflows/fetchData/main.nf"
include {alignGenomes} from "./subworkflows/alignData/main.nf"
include {runScreen;runSNPPipeline;runLocusPipeline} from "./subworkflows/snpdiffs/main.nf"
include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    if(params.runmode == "conda_init"){
        conda_init()
    } else{
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
        if((params.runmode == "align") || (params.runmode == "screen") || (params.runmode == "snp") || (params.runmode == "locus")){

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
                } else if(params.runmode == "locus"){
                    to_align = Channel.empty() 
                }
            }

            // If references are provided, align all queries against all references
            else{
                if ((params.fasta == "") && (params.reads == "")){
                    to_align = Channel.empty()
                } else{
                    to_align = query_data
                    .combine(reference_data)
                    .filter{it -> (it[1].toString() != "null") && (it[3].toString() != "null")} // Can't align without FASTA
                }
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
            } else if(params.runmode == "locus"){
                runLocusPipeline(all_snpdiffs)
            }
        }
    }
}

// Dummy process to stimulate conda env generation
process conda_init {
    executor = 'local'
    cpus = 1
    maxForks = 1
    
    script:
    """
    """
}