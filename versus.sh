#!/bin/sh
TIMESTAMP=`date "+%Y%m%d"`

pwd=`pwd`
mkdir -p $pwd/data/proteomes

### Prepare input files
# Download the ClinVarVariationRelease.xml.gz
cd $pwd/data
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz
gunzip ./ClinVarVariationRelease_00-latest.xml.gz

# Download the fasta sequences
cd $pwd/data/proteomes
for i in {1..100}
do
   wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.$i.protein.faa.gz
done

# Download the tsv file of a list of human enzymes from Uniprot
# cd $pwd/data
# python getHumanEnzymesFromUniprot.py

### Execute the VersUS
cd $pwd/src
python VersUS.py -i config.conf -n $TIMESTAMP

### Register the output to the MySQL database
# pythoh versusDBsql.py ...
