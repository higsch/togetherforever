[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

# Precision Proteogenomics Pipeline
Takes `gtf` files, outputs peptide expression profiles.
Built in cooperation with [Husen's](https://github.com/husensofteng) 
[proteogenomics code](https://github.com/husensofteng/ProteoGenomics).

## Run command for testing purposes
Assumes a lot of dependencies, but docker image will come soon.
```
nextflow main.nf --gtfs "data/gtf/*.gtf" --genome resources/Homo_sapiens.GRCh37.dna.primary_assembly.fa --canonical resources/Homo_sapiens.GRCh38.ENS90.pep.all.fa --mzmldef data/mzml/mzmldef.txt --normpsms data/normpsms/test.txt -resume
```
