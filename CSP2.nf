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

// Set headers
assembly_header = "Isolate_ID\tRead_Type\tRead_Location\tAssembly_Path\n"
isolate_data_header = "Isolate_ID\tIsolate_Type\tAssembly_Path\tContig_Count\tAssembly_Bases\tN50\tN90\tL50\tL90\tSHA256\n"
screen_header = "Query_ID\tReference_ID\tScreening_Category\tCSP2_Screen_SNPs\tQuery_Percent_Aligned\tReference_Percent_Aligned\tQuery_Contig_Count\tQuery_Assembly_Bases\tReferenence_Contig_Count\tReference_Assembly_Bases\tRaw_SNPs\tPurged_Length\tPurged_Identity\tPurged_LengthIdentity\tPurged_Invalid\tPurged_Indel\tPurged_Duplicate\tPurged_Het\tPurged_Density\tFiltered_Query_Edge\tFiltered_Ref_Edge\tFiltered_Both_Edge\tKmer_Similarity\tShared_Kmers\tQuery_Unique_Kmers\tReference_Unique_Kmers\tMUMmer_gSNPs\tMUMmer_gIndels\n"

// Set MUMmer and SNP directories
mummer_directory = file("${output_directory}/MUMmer_Output")
snpdiffs_directory = file("${output_directory}/snpdiffs")
snp_directory = file("${output_directory}/SNP_Analysis")

// Set paths for output files 
snpdiffs_summary_file = file("${output_directory}/Raw_Alignment_Summary.tsv")
isolate_data_file = file("${output_directory}/Isolate_Data.tsv")
screening_results_file = file("${output_directory}/Screening_Results.tsv")

// In --runmode assembly, results save to output_directory
if(params.runmode == "assemble"){
    ref_mode = false

    log_directory = file("${output_directory}")
    assembly_directory = file("${output_directory}")
    assembly_log = file("${log_directory}/Assembly_Data.tsv")

    // Set paths for log files
    user_snpdiffs_list = file("${log_directory}/Imported_SNPDiffs.txt")
    snpdiffs_list_file = file("${log_directory}/All_SNPDiffs.txt")
    screen_log_dir = file("${log_directory}/Screening_Logs")
    snp_log_dir = file("${log_directory}/SNP_Logs")
    
    assembly_log.write("${assembly_header}")
} else{

    // Get reference mode
    if(params.ref_reads == "" && params.ref_fasta == "" && params.ref_id == ""){
        ref_mode = false
    } else{
        ref_mode = true
    }
    
    log_directory = file("${output_directory}/logs")
    assembly_directory = file("${output_directory}/Assemblies")
    assembly_log = file("${log_directory}/Assembly_Data.tsv")

    // Set paths for log files
    user_snpdiffs_list = file("${log_directory}/Imported_SNPDiffs.txt")
    snpdiffs_list_file = file("${log_directory}/All_SNPDiffs.txt")
    screen_log_dir = file("${log_directory}/Screening_Logs")
    snp_log_dir = file("${log_directory}/SNP_Logs")

    log_directory.mkdirs()
    mummer_directory.mkdirs()
    snpdiffs_directory.mkdirs()

    // Establish Isolate_Data.tsv
    isolate_data_file.write("${isolate_data_header}")

    // If --reads/--ref_reads are provided, prepare a directory for assemblies
    if((params.reads != "") || (params.ref_reads != "")){
        assembly_directory.mkdirs()
        assembly_log.write("${assembly_header}")    
    }

    // If runmode is snp, prepare a directory for SNP analysis + logs
    if(params.runmode == "snp"){
        snp_directory.mkdirs()
        snp_log_dir.mkdirs()
    }

    // If runmode is screen, prepare a directory for screening logs
    if(params.runmode == "screen"){
        screen_log_dir.mkdirs()
        
        // Establish Screening_Results.tsv
        screening_results_file.write("${screen_header}")
    }
}

// Parameterize variables to pass between scripts
params.output_directory = file(output_directory)
params.log_directory = file(log_directory)
params.assembly_directory = file(assembly_directory)

params.assembly_log = file(assembly_log)

params.mummer_directory = file(mummer_directory)

params.snpdiffs_directory = file(snpdiffs_directory)
params.snpdiffs_list_file = file(snpdiffs_list_file)
params.snpdiffs_summary_file = file(snpdiffs_summary_file)
params.user_snpdiffs_list = file(user_snpdiffs_list)

params.snp_directory = file(snp_directory)
params.isolate_data_file = file(isolate_data_file)
params.screening_results_file = file(screening_results_file)

params.screen_log_dir = file(screen_log_dir)
params.snp_log_dir = file(snp_log_dir)

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
include {runScreen} from "./subworkflows/snpdiffs/main.nf"

//include {saveIsolateLog} from "./subworkflows/logging/main.nf"
//include {runScreen;runSNPPipeline} from "./subworkflows/snpdiffs/main.nf"
//include {runRefChooser} from "./subworkflows/refchooser/main.nf"

workflow{
    
    // Read in data
    input_data = fetchData()

    // Create channel for pre-aligned data
    already_aligned = input_data.snpdiff_data
    .map { it -> tuple([it[0], it[1]].sort().join(','), it[2]) }
    .collect().flatten().collate(2)
    
    // If run mode is 'assemble', tasks are complete
    if((params.runmode == "align") || (params.runmode == "screen")){

        // If there is no reference data, align all query_data against each other
        if(!ref_mode){
            
            seen_combinations = []
            
            to_align = input_data.query_data.combine(input_data.query_data) // Self-combine query data
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
        
        // If references are provided, align all queries against all references
        else{
            to_align = input_data.query_data
            .combine(input_data.reference_data)
            .filter{it -> (it[1].toString() != "null") && (it[3].toString() != "null")} // Can't align without FASTA
        }

        // Don't align things that are already aligned via --snpdiffs
        mummer_results = to_align
        .map { it -> tuple([it[0], it[2]].sort().join(','),it[0], it[1], it[2], it[3]) }
        .join(already_aligned,by:0,remainder:true)
        .filter{it -> it[5].toString() == "null"}
        .map{it -> [it[1], it[2], it[3], it[4]]}
        | alignGenomes | collect | flatten | collate(3)
        
        all_snpdiffs = input_data.snpdiff_data.concat(mummer_results)
        .unique{it->it[2]}
        .collect().flatten().collate(3)
        .ifEmpty { error "No .snpdiffs to process..." }
        
        // Save raw alignment log
        all_snpdiffs.collect{it[2]} | saveMUMmerLog 
        
        if(params.runmode == "screen"){
                        
            // If references are provided, ensure that all snpdiffs are processed in Query-Reference order
            // If no references are provided, process snpdiffs as is
            if(ref_mode){

                ref_list = input_data.reference_data
                    .unique{it -> it[0]}
                    .map{it -> "${it[0]}"}
                    .collect()

                query_list = input_data.query_data
                    .unique{it -> it[0]}
                    .map{it -> "${it[0]}"}
                    .collect()

                all_snpdiffs = all_snpdiffs
                .branch{it ->
                    forward: query_list.contains("${it[0]}") && ref_list.contains("${it[1]}")
                        return(tuple([it[0], it[1]].sort().join(','),it[0],it[1],it[2]))
                    reverse: ref_list.contains("${it[0]}") && query_list.contains("${it[1]}")
                        return(tuple([it[0], it[1]].sort().join(','),it[1],it[0],it[2]))
                    other: true
                        return(tuple(it[0],it[1],it[2]))}

                // Drop snpdiffs that are in the reverse orientation if a forward orientation exists
                retained_reverse = all_snpdiffs.reverse
                .join(all_snpdiffs.forward, by:0, remainder:true)
                .filter{it -> it[4].toString() == "null"}
                .map{it -> [it[1], it[2], it[3]]}

                
                all_snpdiffs.forward
                .map{it -> [it[1], it[2], it[3]]}
                .concat(retained_reverse)
                .concat(all_snpdiffs.other)
                .unique{it -> it[2]}.collect().flatten().collate(3)
                | runScreen

                } else{
                    all_snpdiffs.collect().flatten().collate(3)
                    | runScreen
                }
        }
    } else if(params.runmode == "snp"){
        print("SNP")
    }
}

process saveMUMmerLog{
// Takes: Flattened list of snpdiffs files and generates a TSV report via python
// Returns: Path to list of snpdiffs files
    executor = 'local'
    cpus = 1
    maxForks = 1

    input:
    val(snpdiffs_paths)

    script:
    saveSNPDiffs = file("$projectDir/bin/saveSNPDiffs.py")
    snpdiffs_list_file.write(snpdiffs_paths.join('\n') + '\n')
    """
    $params.load_python_module
    python $saveSNPDiffs "${snpdiffs_list_file}" "${snpdiffs_summary_file}"
    """
}

        
        /*
        else if(params.runmode == "screen"){

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
    else{
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

            if(params.runmode == "screen"){
                input_data.snpdiffs_data.concat(mummer_results) | runScreen // Compare snpdiffs to generate a summary
            }
            else if(params.runmode == "snp"){
                input_data.snpdiffs_data.concat(mummer_results) | runSNPPipeline // Generate pairwise SNP distances and alignments against each reference
            }
        } 
    }
}








*/