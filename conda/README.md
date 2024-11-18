To use the YAML file:

```
conda env create -f environment.yaml
conda activate csp2_env
nextflow run CSP2.nf -profile local_multithread --cores XX --runmode XX ...
```
