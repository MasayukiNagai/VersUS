# test script
from main import *
from XMLHandler import *

genes_file = './data/gene/HumanEnzWithEC.csv'
genes_dict = readHumanGenesEC(genes_file)

clinvar_file = './data/clinvar/ClinVarVariationRelease_00-latest_weekly.xml'
clinvar_file = './data/clinvar/clinvar_tail_5.xml'
out_file = './data/VUS.csv'
print('---------- VariationParser ----------')
data = readClinVarVariationsXML(clinvar_file, out_file, genes_dict)


