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
