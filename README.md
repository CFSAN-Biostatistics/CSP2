# Yenta (v.0.5)  

Yenta is a Nextflow pipeline for fast and accurate SNP distance estimation from WGS read data or genome assemblies.  

---
## Software Dependencies  
The following software are required to run Yenta. Software version used during Yenta development noted in parentheses.  

- Nextflow (22.10.7)  
- Python (3.8.1)  
- BEDTools (2.26.0)  
- MUmmer (4.0.0)  
- SKESA (2.5.0) [Only required if starting from raw reads]  
  
---
## Installation  
Yenta can be run by cloning the GitHub repo.  

```
git clone https://github.com/CFSAN-Biostatistics/Yenta.git
```

---
## Configuration  
Yenta options can be specified on the command line, or through the Nextflow configuration files [nexflow.config + conf/profiles.config]  

**Options with defaults include**:  
| Parameter    | Description                                                                                                               | Default                                   |
|--------------|---------------------------------------------------------------------------------------------------------------------------|-------------------------------------------|
| cores        | CPUs per node                                                                                                             | 1                                         |
| align_cov    | Only consider queries where either the query or reference genome are covered by at least <align_cov>%                     | 85                                        |
| ref_iden     | Only consider SNPs from contig alignments with <ref_iden>% identity                                                       | 99                                        |
| ref_edge     | Remove SNPs that occur within <ref_edge>bp from the end of the reference contig                                           | 500                                       |
| query_edge   | Remove SNPs that occur within <query_edge>bp from the end of the query contig                                             | 500                                       |
| out          | Path to output folder (Must not exist)                                                                                    | ./YENTA_${new java.util.Date().getTime()} |
| forward      | Suffix for forward query reads                                                                                            | _1.fastq.gz                               |
| ref_forward  | Suffix for forward reference reads                                                                                        | _1.fastq.gz                               |
| reverse      | Suffix for reverse query reads                                                                                            | _2.fastq.gz                               |
| ref_reverse  | Suffix for reverse reference reads                                                                                        | _2.fastq.gz                               |
| readext      | Extension for query reads                                                                                                 | fastq.gz                                  |
| ref_readext  | Extension for reference reads                                                                                             | fastq.gz                                  |
| readtype     | Query read naming convention. Options include srazip (_1/2.fastq.gz), illumina (_R1/2_001.fastq.gz), sra (_1/2.fastq)     | srazip                                    |
| ref_readtype | Reference read naming convention. Options include srazip (_1/2.fastq.gz), illumina (_R1/2_001.fastq.gz), sra (_1/2.fastq) | srazip                                    |


**Options without defaults include**:  
| Parameter       | Description                                                                                                       | Default |
|-----------------|-------------------------------------------------------------------------------------------------------------------|---------|
| reads           | Location of query read data (Path to directory, or path to file with multiple directories)                        | USER    |
| fasta           | Location of query assembly data (Path to directory containing FASTAs, path to FASTA, path to multiple FASTAs)     | USER    |
| ref_reads       | Location of reference read data (Path to directory, or path to file with multiple directories)                    | USER    |
| ref_fasta       | Location of reference assembly data (Path to directory containing FASTAs, path to FASTA, path to multiple FASTAs) | USER    |
| python_module   | Name of Python module if 'module load PYTHON' statement is required.                                                      | USER    |
| mummer_module   | Name of MUmmer module if 'module load MUMMER' statement is required.                                                      | USER    |
| skesa_module    | Name of SKESA module if 'module load SKESA' statement is required.                                                        | USER    |
| bedtools_module | Name of BEDTools module if 'module load BEDTOOLS' statement is required.                                                  | USER    |

---

## Example Runs
Here are a few examples of how you can run Yenta:  

```
module load nextflow
git clone https://github.com/CFSAN-Biostatistics/Yenta.git
cd Yenta

# Assemble reads for queries and compare to single reference fasta using 20 CPUs, output to ./Tiny_Test
nextflow run Yenta.nf                           // Run Yenta  
--reads assets/test_reads                       // Use all reads from this folder  
--readtype srazip                               // Reads are zipped SRA (_1/2.fastq.gz)  
--ref_fasta assets/test_ref/SRR10831135.fasta   // Compare query to single reference  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test                                 // Save to ./Tiny_Test  

# Compare assembled queries to multiple references fasta using 20 CPUs, output to ./Tiny_Test_2, use specific python module
nextflow run Yenta.nf                           // Run Yenta  
--fasta assets/test_fasta                       // Use all assemblies from this folder  
--ref_fasta assets/test_ref                     // Compare queries to this assembled reference and...  
--ref_reads assets/test_ref                     // Assemble and use this reference too  
--ref_readtype srazip                           // Reads are zipped SRA (_1/2.fastq.gz)  
--cores 20                                      // Use 20 CPUs per node  
--out Tiny_Test_2                               // Save to ./Tiny_Test_2  
--python_module python/3.8.1                    // run 'module load python/3.8.1' in scripts using Python 

# Adjust filtering criteria  
nextflow run Yenta.nf                           // Run Yenta  
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
Yenta outputs basic information about reference and query isolates (contig counts, base counts, data location), as well as detailed alignment and SNP information about each query/reference pair. For example, this is the output from the top command above.  

### /your/out/dir/Sample_Data.tsv
| Sample_ID  | Data_Type | Read_Data                                                                                                                                                                                                             | Assembly_Data                                                                                                        | Assembly_Contigs | Assembly_Bases |
|------------|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|------------------|----------------|
| SRR1849356 | Paired    | assets/test_reads/SRR1849356_1.fastq.gz;assets/test_reads/SRR1849356_2.fastq.gz | SRR1849356/SRR1849356.fasta | 63               | 2910691        |

### /your/out/dir/Reference_Strain_Data/Reference_Data.tsv
| Reference_ID | Data_Type | Read_Data | Assembly_Data                                                                                        | Assembly_Contigs | Assembly_Bases |
|--------------|-----------|-----------|------------------------------------------------------------------------------------------------------|------------------|----------------|
| SRR10831135  | Assembly  | NA        | assets/test_ref/SRR10831135.fasta | 47               | 2914526        |

### /your/out/dir/Screening_Results/MUmmer_DNADiff_Results.tsv
| Sample_ID  | Reference_ID | Percent_Reference_Covered | Percent_Query_Covered | Category | Yenta_SNPs | gSNPs | Filtered_Identity | Filtered_Edge | Filtered_Duplicated | Rejected_Density_1000 | Rejected_Density_125 | Rejected_Density_15 |
|------------|--------------|---------------------------|-----------------------|----------|------------|-------|-------------------|---------------|---------------------|-----------------------|----------------------|---------------------|
| SRR1849356 | SRR10831135  | 99.78137096735455         | 99.91287292261528     | PASS     | 3          | 1     | 1                 | 0             | 0                   | 0                     | 0                    | 0                   |

