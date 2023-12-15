# CFSAN SNP Pipeline v2 (CSP2)  

**Important Note:** *CSP2 is currently under development, and has not been validated for non-research purposes. Current workflows and data processing parameters may change prior to full release version.*

### CSP2 is a Nextflow pipeline for fast and accurate microbial SNP distance estimation

CSP2 can be run in one of two main modes:  

#### 1) "Screening Mode" (--runmode screen)
Given a set of reference isolates provided by the user, screen query assemblies and generate a summary report.  
   
   *e.g.* Screen incoming NGS results against a database of known lab controls to quickly identify potential contamination. 

   ```
   nextflow run CSP2.nf --runmode screen --ref_fasta /path/to/lab/control/assemblies --reads /path/to/sequencing/run/fastq --out contam_check
   ```

   This command will:
   1) Detect read files
   2) Generate assemblies using SKESA
   3) Screen each query assembly against each lab control sequence provided
   4) Provide a report with alignment overlap statistics and SNP distances for each query/reference combination. 

#### 2) "SNP Pipeline Mode" (--runmode snp)
Generate the pairwise SNP distances and alignments for a cluster of isolates using either user-provided references or reference isolates automatically identified by RefChooser. 

**Note:** Testing is underway to determine how the underlying cluster diversity impacts distances estimates. Current comparisons are based on SNP clusters with densities in the 0 ~ 100 SNP range.  

   *e.g.* Get the pairwise distances for a group of isolates all found in the same soil sample

   ```
   nextflow run CSP2.nf --runmode snp --fasta /path/to/soil/isolate/assemblies
   ```
   This command will:
   1) Detect query assemblies
   2) Choose a reference isolate via RefChooser
   3) Align each query isolate to the reference
   4) Summarize individual mapping results and generate a full pairwise SNP matrix, SNP alignment, and SNP coordinate data

OR 
   ```
   nextflow run CSP2.nf --runmode snp --fasta /path/to/soil/isolate/assemblies --ref_fasta /path/to/favorite/soil/assembly.fa
   ```
   This command will:
   1) Detect query assemblies
   2) Align each query isolate to /path/to/favorite/soil/assembly.fa
   3) Summarize individual mapping results and generate a full pairwise SNP matrix, SNP alignment, and SNP coordinate data
   
---
## Software Dependencies  
The following software are required to run CSP2. Software version used during CSP2 development noted in parentheses.  

- Nextflow (22.10.7)  
- Python (3.8.1)  
  - pybedtools  
- BEDTools (2.26.0)  
- MUmmer (4.0.0)  
- SKESA (2.5.0) [Only required if starting from raw reads]  
  
---
## Installing CSP2 
CSP2 can be installed by cloning the GitHub repo.  

```
git clone https://github.com/CFSAN-Biostatistics/CSP2.git
```
---
## Tips for configuring CSP2  
CSP2 options can be specified on the command line, or through the Nextflow configuration files detailed in the next section. Feel free to skip this section if you're familiar with editing Nextflow configuration files.  

There are two main configuration files associated with CSP2:  
- [nextflow.config](nextflow.config)  
- [profiles.config](conf/profiles.config)
- **Note:** Any options set in either of these files are overruled by options set on the command line  

The [nextflow.config](nextflow.config) file shown below can be edited to include 'hard-set' parameters you want to be used for every CSP2 run.   

```
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CSP2 Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Default config options for all compute environments
----------------------------------------------------------------------------------------
*/

// Import profile settings
includeConfig "${projectDir}/conf/profiles.config"

// Global default params
params {

    // Setting output directory 

    // Set name for output folder/file prefixes
    out = "CSP2_${new java.util.Date().getTime()}"

    // Set output parent directory [Default: CWD; Set this to have all output go to the same parent folder, with unique IDs set by --out]
    outroot = ""

    
    // CSP2 can run in the following run-modes:

    // assemble: Assemble read data (--reads/--ref_reads) into FASTA via SKESA (ignores --fasta/--ref_fasta/--snpdiffs)
    // align: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta), run MUMmer alignment analysis for each query/ref combination (ignores --snpdiffs)
    // screen: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta) and/or MUMmer output (.snpdiffs), create a report for raw SNP distances between each query and reference assembly
    // snp: Given query data (--reads/--fasta) and reference data (--ref_reads/--ref_fasta) and/or MUMmer output (.snpdiffs), generate alignments and pairwise distances for all queries based on each reference dataset
    
    runmode = ""

    // Location for isolate sequence data
    reads = ""
    fasta = ""

    // Location for reference sequence data
    ref_reads = ""
    ref_fasta = ""

    // Location for snpdiffs files
    snpdiffs = ""
    
    // Read read_info
    readext = "fastq.gz"
    forward = "_1.fastq.gz"
    reverse = "_2.fastq.gz"

    ref_readext = "fastq.gz"
    ref_forward = "_1.fastq.gz"
    ref_reverse = "_2.fastq.gz"

    // Analytical variables

    // Only consider queries if the reference genome is covered by at least <min_cov>% [Default: 85]
    min_cov = 85

    // Only consider SNPs from contig alignments longer than <min_len> bp [Default: 500]
    min_len = 500

    // Only consider SNPs from contig alignments with <min_iden>% identity [Default: 99]
    min_iden = 99

    // Remove SNPs that occur within <ref_edge>bp from the end of the reference contig [Default: 250]
    ref_edge = 250

    // Remove SNPs that occur within <query_edge>bp from the end of the query contig [Default: 250]
    query_edge = 250

    // SNP density filters: Given density windows provided by dwin, purge windows where more than the allowable window SNPs (wsnps) are found
    // Default: 3 max per 1000bp, 2 max per 125bp, 1 max per 15bp, filtered from biggest window to smallest
    dwin = "1000,125,15"
    wsnps = "3,2,1"

    // If running refchooser in snp mode, compare queries to the top X references [Default: 1]
    n_ref = 1

    // If the assembly file contains the string <trim_name>, remove it from the sample name (e.g. '_contigs_skesa')
    trim_name = '""'
}
```

The [profiles.config](conf/profiles.config) file shown below can be edited to add information about your computing environment. These include:  

- Number of cores per cpu (cores)  
- Number of concurrent processess allowed at once (process.maxForks)  
- Names of modules to load for Python, MUmmer, SKESA, and BEDTools (if required)  
  - If modules are specified, they are loaded automatically by Nextflow via  ```module load -s <MODULE>```

An example configuration setup (slurmHPC) is provided as a model.  

```
profiles {
    standard {
        process.executor = 'local'
        process.cpus = 1
        process.maxForks = 1
        params.python_module = ""
        params.mummer_module = ""
        params.skesa_module = ""
        params.bedtools_module = ""
    }
    local_multithread {
        process.executor = 'local'
        process.maxForks = 1
        params.cores = 1
        process.cpus = "${params.cores}"
        params.python_module = ""
        params.mummer_module = ""
        params.skesa_module = ""
        params.bedtools_module = ""
    }
    slurmHPC {
        process.executor = 'slurm'
        process.maxForks = 30
        params.cores = 20
        process.cpus = "${params.cores}"
        params.python_module = "python/3.8.1"
        params.mummer_module = "mummer/4.0.0"
        params.skesa_module = "skesa/2.5.0"
        params.bedtools_module = "bedtools"
    }
}
```
---


**Options with defaults include**:  
| Parameter     | Description                                                                                                | Default Value                          |
|---------------|------------------------------------------------------------------------------------------------------------|----------------------------------------|
| --outroot         | Base directory to create output folder                                                       | $CWD |
| --out         | Name of the output folder to create (must not exist)                                                       | CSP2_${new java.util.Date().getTime()} |
| --forward     | Full file extension for forward/left reads of query                                                        | _1.fastq.gz                            |
| --reverse     | Full file extension for reverse/right reads of reference                                                   | _2.fastq.gz                            |
| --ref_forward | Full file extension for forward/left reads of reference                                                    | _1.fastq.gz                            |
| --ref_reverse | Full file extension for reverse/right reads of reference                                                       | _2.fastq.gz                            |
| --readext     | Extension for single-end reads for query                                                                   | fastq.gz                               |
| --ref_readext | Extension for single-end reads for reference                                                               | fastq.gz                               |
| --min_cov     | Do not analyze queries that cover less than <min_cov>% of the reference assembly                           | 85                                     |
| --min_iden    | Only consider alignments where the percent identity is at least <min_iden>%                                | 99                                     |
| --min_len     | Only consider alignments that span at least <min_len>bp                                                    | 500                                    |
| --dwin        | A comma-separated list of windows to check SNP densities                                                   | 1000,125,15                            |
| --wsnps       | The maximum number of SNPs allowed in the corresponding window from --dwin                                 | 3,2,1                                  |
| --query_edge  | Only consider SNPs that occur within <query_edge>bp of the end of a query contig                           | 250                                    |
| --ref_edge    | Only consider SNPs that occur within <query_edge>bp of the end of a reference contig                       | 250                                    |
| --n_ref       | The number of RefChooser reference isolates to consider (only applied if using RefChooser)                 | 1                                      |

**Options without defaults include**:  
| Parameter           | Description                                                                                                       |
|---------------------|-------------------------------------------------------------------------------------------------------------------|
| --reads             | Location of query read data (Path to directory, or path to file with multiple directories)                        |
| --fasta             | Location of query assembly data (Path to directory containing FASTAs, path to FASTA, path to multiple FASTAs)     |
| --ref_reads         | Location of reference read data (Path to directory, or path to file with multiple directories)                    |
| --ref_fasta         | Location of reference assembly data (Path to directory containing FASTAs, path to FASTA, path to multiple FASTAs) |
| --python_module     | Name of Python module if 'module load PYTHON' statement is required.                                              |
| --mummer_module     | Name of MUmmer module if 'module load MUMMER' statement is required.                                              |
| --skesa_module      | Name of SKESA module if 'module load SKESA' statement is required.                                                |
| --refchooser_module | Name of RefChooser module if 'module load REFCHOOSER' statement is required.                                      |
| --bedtools_module   | Name of BEDTools module if 'module load BEDTOOLS' statement is required.                                          |
| --trim_name         | A string in assembly file names that you want to remove from sample IDs (e.g., _contigs_skesa)                    |

---


## Examples

Here are a few examples of how you can use CSP2: 

1. You want to know whether incoming read data is actually from the lab control strain. 

```


# Assemble reads for queries and compare to single reference fasta using 20 CPUs, output to ./Tiny_Test
nextflow run /path/to/CSP2.nf                   // Run CSP2  
--reads assets/test_reads                       // Use all reads from this folder  
--readtype srazip                               // Reads are zipped SRA (_1/2.fastq.gz)  
--ref_fasta assets/test_ref/SRR10831135.fasta   // Compare query to single reference  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test                                 // Save to ./Tiny_Test  

# Compare assembled queries to multiple references fasta using 20 CPUs, output to ./Tiny_Test_2, use specific python module
nextflow run CSP2.nf                           // Run CSP2  
--fasta assets/test_fasta                       // Use all assemblies from this folder  
--ref_fasta assets/test_ref                     // Compare queries to this assembled reference and...  
--ref_reads assets/test_ref                     // Assemble and use this reference too  
--ref_readtype srazip                           // Reads are zipped SRA (_1/2.fastq.gz)  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test_2                               // Save to ./Tiny_Test_2  
--python_module python/3.8.1                    // run 'module load python/3.8.1' in scripts using Python 

# Adjust filtering criteria  
nextflow run CSP2.nf                           // Run CSP2  
--reads assets/test_reads                       // Use all reads from this folder  
--readtype srazip                               // Reads are zipped SRA (_1/2.fastq.gz)  
--ref_fasta assets/test_ref/SRR10831135.fasta   // Compare query to single reference  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test_3                               // Save to ./Tiny_Test_3  
--align_cov 95                                  // Only consider query/ref pairs where 95% of each genome is aligned   
--ref_edge 250                                  // Remove SNPs that are 250bp from the edge of a reference contig  
--query_edge 600                                // Remove SNPs that are 600bp from the edge of a query contig  
--ref_iden 95                                   // Allow SNPs from contig alignments with at least 95% sequence identity  
```

---
## Output
CSP2 outputs basic information about reference and query isolates (contig counts, base counts, data location), as well as detailed alignment and SNP information about each query/reference pair. For example, this is the output from the top command above.  

### /your/out/dir/Sample_Data.tsv
| Sample_ID  | Data_Type | Read_Data                                                                                                                                                                                                             | Assembly_Data                                                                                                        | Assembly_Contigs | Assembly_Bases |
|------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|------------------|----------------|
| SRR1849356 | Paired    | assets/test_reads/SRR1849356_1.fastq.gz;assets/test_reads/SRR1849356_2.fastq.gz | SRR1849356/SRR1849356.fasta | 63               | 2910691        |

### /your/out/dir/Reference_Strain_Data/Reference_Data.tsv
| Reference_ID | Data_Type | Read_Data | Assembly_Data                                                                                        | Assembly_Contigs | Assembly_Bases |
|--------------|-----------|-----------|------------------------------------------------------------------------------------------------------|------------------|----------------|
| SRR10831135  | Assembly  | NA        | assets/test_ref/SRR10831135.fasta | 47               | 2914526        |

### /your/out/dir/Screening_Results/MUmmer_DNADiff_Results.tsv
| Sample_ID  | Reference_ID | Percent_Reference_Covered | Percent_Query_Covered | Category | CSP2_SNPs | gSNPs | Filtered_Identity | Filtered_Edge | Filtered_Duplicated | Rejected_Density_1000 | Rejected_Density_125 | Rejected_Density_15 |
|------------|--------------|---------------------------|-----------------------|----------|------------|-------|-------------------|---------------|---------------------|-----------------------|----------------------|---------------------|
| SRR1849356 | SRR10831135  | 99.78137096735455         | 99.91287292261528     | PASS     | 3          | 1     | 1                 | 0             | 0                   | 0                     | 0                    | 0                   |

