# test script
from main import *
from XMLHandler import *
from Processor import * 

processor = Processor()

print('---------- read Gene File ----------')
genes_file = './data/gene/HumanEnzWithEC.csv'
genes_dict = readHumanGenesEC(genes_file)

# print('---------- VariationParser ----------')
# clinvar_file = './data/clinvar/ClinVarVariationRelease_00-latest_weekly.xml'
# # clinvar_file = './data/clinvar/clinvar_tail_5.xml'
# vus_dict = readClinVarVariationsXML(clinvar_file, genes_dict)

# print('---------- write to CSV ----------')
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'pos')
outfile_path = './data/VUS.csv'
# processor.write_to_csv(vus_dict, header, outfile_path)

print('---------- read CSV ----------')
csv_path = outfile_path
vus_dict = processor.read_csv_to_dict(csv_path)
print(vus_dict[0])

print('---------- make sequence dict ----------')
seq_dir_path = './data/fasta_sequences'
seq_dict = processor.make_seq_dict(seq_dir_path)

print('---------- add seq to vus_dict ----------')
seq_range = 12
vus_dict, unfound_seq = processor.add_seq_to_dict(vus_dict, seq_dict, seq_range)
print(unfound_seq)