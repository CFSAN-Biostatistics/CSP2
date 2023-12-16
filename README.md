<p align="center">
  <img src="CSP2_Logo_B.jpg" alt="drawing" width="250"/>
</p>

# CFSAN SNP Pipeline 2 (CSP2)  


**Important Note:** *CSP2 is currently under development, and has not been validated for non-research purposes. Current workflows and data processing parameters may change prior to full release version.*

## CSP2 is a Nextflow pipeline for rapid, accurate SNP distance estimation from assembly data  

All CSP2 sequence comparisons happen at the assembly level, but if reads are provided CSP2 will perform a genome assembly using *SKESA*. CSP2 has two main run modes:  

#### 1) "Screening Mode" (*--runmode screen*)
Given one or more user-provided reference isolates (*--ref_reads*; *--ref_fasta*), get alignment statistics and SNP distances for query isolates (*--reads*; *--fasta*)
   
#### 2) "SNP Pipeline Mode" (*--runmode snp*)
Generate pairwise SNP distances and alignments for 2+ isolates (*--reads*; *--fasta*) based on comparisons to:  
- One or more user-provided references (*--ref_reads*; *--ref_fasta*), or  
- One or more reference isolates selected by RefChooser (*--n_ref*)

**Important Note**: *Testing is underway to determine how the underlying cluster diversity impacts distances estimates. Current comparisons are based on strains clusters with SNP densities in the 0 ~ 150 SNP range.*

In either case, CSP2 calls MUMmer for alignment and if a sufficient portion of the reference genome is aligned (*--min_cov*), that data is passed through a set of filters, including the automated removal of:  
- Sites from very short (*--min_len*) or poorly aligned contigs (*--min_iden*)
- Multiply aligned sites
- Sites from regions of high SNP density (*--dwin*/*--wsnps*)
- Non-base sites (e.g., 'N' or '?')
- Heterozygous sites
- Sites close to the contig edge (*--query_edge*/*--ref_edge*)
- Indels (**for now**)

This final dataset is summarized into a *.snpdiffs* file, which contains:  
1. A one-line header with alignment statistics  
2. A BED file of contig mappings that pass QC  
3. Information about SNPs (if present)

To avoid unnecessary realignment, once a .snpdiffs file is generated under a particular set of QC parameters (which is hardcoded into the file) these files can be used in other CSP2 runs via the *--snpdiffs* argument (if using the same QC parameters). 

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
CSP2 can be installed by cloning the GitHub repo and configuring the [nextflow.config](nextflow.config) and [profiles.config](conf/profiles.config) to suit your needs  

```
git clone https://github.com/CFSAN-Biostatistics/CSP2.git
```

## Tips for configuring CSP2  
CSP2 options can be specified on the command line, or through the Nextflow configuration files detailed in the next section. Feel free to skip this section if you're familiar with editing Nextflow configuration files.  

There are two main configuration files associated with CSP2:  

- The profiles.config file is where you add custom information about your computing environment. An example configuration setup (slurmHPC) is provided as a model.

```
profiles {
    standard {
        process.executor = 'local'
        params.cores = 1
        params.python_module = ""
        params.mummer_module = ""
        params.skesa_module = ""
        params.bedtools_module = ""
        params.refchooser_module = ""
    }
    slurmHPC {
        process.executor = 'slurm'
        params.cores = 20
        params.python_module = "python/3.8.1"
        params.mummer_module = "mummer/4.0.0"
        params.skesa_module = "skesa/2.5.0"
        params.bedtools_module = "bedtools"
        params.refchooser_module = "refchooser/0.2.1"
    }
}
```
- If you add your own profile, be sure to note it on the command line (one hypen)
```
nextflow run CSP2.nf -profile myNewProfile <args>
```
- To save an HTML report from the Nextflow run:
```
nextflow run CSP2.nf -with-report myReport.html <args>
```


- The nextflow.config file is where you can change other aspects of the CSP2 run, including data location, QC parameters, and all the options listed below:

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

The repo contains small test datasets to ensure things are running as expected. Here are a few examples of how you can use CSP2 in screening mode or in SNP pipeline mode. 

**Screening Mode**

**SNP Pipeline Mode**


```
nextflow run /path/to/CSP2.nf                   // Run CSP2  
-profile myHPC                                  // Choose run profile
-with-report Contamination.html                 // Save an HTML report
--runmode screen                                // Compare each query to the reference (Don't need pairwise between queries)
--forward _1.fq.gz                              // Forward reads don't match the default '_1.fastq.gz'
--reverse _2.fq.gz                              // Reverse reads don't match the default '_2.fastq.gz'
--readext fq.gz                                 // Single-end reads don't match the default 'fastq.gz'
--ref_fasta assets/test_ref/SRR10831135.fasta   // Compare query to single reference  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test                                 // Save to ./Tiny_Test  
```

---
## Output
