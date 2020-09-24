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
# outfile_path = './data/VUS.tsv'
# processor.write_to_tsv(vus_dict, header, outfile_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS.tsv'
# vus_dict = processor.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make sequence dict ----------')
# seq_dir_path = './data/fasta_sequences'
# seq_dict = processor.make_seq_dict(seq_dir_path)
# print(next(iter(seq_dict.keys())))

# print('---------- add seq to vus_dict ----------')
# seq_range = 12
# vus_dict, unfound_seq = processor.add_seq_to_dict(vus_dict, seq_dict, seq_range)
# print(f'NP numbers missing their seq {len(unfound_seq)}')

# print('---------- write to CSV ----------')
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'sequence')
# outfile_path = './data/VUS.tsv'
# processor.write_to_tsv(vus_dict, header, outfile_path)

print('********** read TSV **********')
tsv_path = './data/VUS.tsv'
vus_dict = processor.read_tsv_to_dict(tsv_path)
print(f'The number of VUS: {len(vus_dict)}')
print(vus_dict[0])

print('---------- make a fasta file for blast ----------')
fasta_for_blast_path = './data/blast/vus_blast.fasta'
processor.make_fasta_for_blast(vus_dict, fasta_for_blast_path)
print('Done')

print('---------- run BLAST ----------')
blast_output_path = './data/blast/blast_results.fasta'
processor.blast_locally(fasta_for_blast_path, blast_output_path)