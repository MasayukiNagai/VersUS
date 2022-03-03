#!/bin/sh
TIMESTAMP=`date "+%Y%m%d"`

PWD=`pwd`
mkdir -p $PWD/data/fasta

### Prepare input files
# Download the ClinVarVariationRelease.xml.gz
cd $PWD/data
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz

# Download the fasta sequences
cd $PWD/data/fasta
for i in {1..100}
do
   wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.$i.protein.faa.gz
done

# Download the tsv file of a list of human enzymes from Uniprot
# cd $PWD/data
# python getHumanEnzymesFromUniprot.py

### Execute the VersUS
python VersUS.py -i $PWD/config.conf -n $TIMESTAMP

### Register the output to the MySQL database
# pythoh versusDBsql.py ...
