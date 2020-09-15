# test script
from main import *
from XMLParser import *

genes_file = './data/HumanEnzWithEC.csv'
genes_dict = readHumanGenesEC(genes_file)

clinvar_file = './data/ClinVarVariationRelease_00-latest_weekly.xml'
out_file = './data/VUS.csv'
readClinVarVariationsXML(clinvar_file, out_file, genes_dict)