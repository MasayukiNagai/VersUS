# VersUS
## Overview

VersUS provides a list of promising leads in the area of rare disease research by extracting variants of uncertain significance (VUS) from the ClinVar database and annotating them with enzyme type, PDB codes, gnomAF, and CADD scores. 

## Requirements

### Tools

- Python 3
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
- [Ensemble Variant Effect Predictor (VEP)](https://asia.ensembl.org/info/docs/tools/vep/index.html)

### Python Dependencies

`pip install` the following

* **lxml**
* **Bio**
* **selenium** (for retriving CADD scores)
* **chromedriver_binary** (for retriving CADD scores)
* **wget** (for retriving CADD scores)

## Usage
```bash
$ python VersUS.py -c config.conf
```

### config.conf

```
[general]
## input files
# a list of human enzymes with corresponding EC numbers
genes=../data/HumanEnzWithEC.csv
# file path for a ClinVar Variation Release
variations=../data/ClinVarVariationRelease_00-latest_weekly.xml 
#direcotry path for proteomes fasta sequences
proteomes=../data/fasta_sequences

## command line tools
# directory path for the BLAST+
blast=../ncbi/blast/bin
# direcotry path for the ensemble VEP 
vep=../ensembl-vep

# True if you want to retrieve CADD scores, otherwise False
cadd=True

## output direcotry
outdir=./result
# directory to store intermediate files
intermediates=./intermediates

verbose=True

[parameters]
fasta_window=12
evalue=10.0
```



## Output

```
result
├── versus_logger.log
└── VersUS.csv
intermediates
├── blast_vus_input.fasta
├── blast_vus_results.xml
├── vep_vus_input.tsv
├── vep_vus_results.tsv
├── vep_vus_results.tsv_summary.html
├── cadd_vus_input.tsv
└── cadd_vus_scores.tsv.gz
```



## Acknowledgement

This project has been done under the supervision of Professor. Daniel Gurnon and Professor. Chad Byers at DePauw University. 
