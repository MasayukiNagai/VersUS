# test script
from main import *
from XMLHandler import *
from Processor import * 
from WebsiteHandler import *


processor = Processor()
web_handler = WebsiteHandler()

# print('---------- read Gene File ----------')
# genes_file = './data/gene/HumanEnzWithEC.csv'
# genes_dict = readHumanGenesEC(genes_file)

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
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window')
# outfile_path = './data/VUS.tsv'
# processor.write_to_tsv(vus_dict, header, outfile_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS.tsv'
# vus_dict = processor.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make a fasta file for blast ----------')
# fasta_for_blast_path = './data/blast/vus_blast.fasta'
# processor.make_fasta_for_blast(vus_dict, fasta_for_blast_path)
# print('Done')

# print('---------- run BLAST ----------')
# blast_output_path = './data/blast/blast_results.xml'
# processor.blast_locally(fasta_for_blast_path, blast_output_path)

# print('---------- read BLAST results ----------')
# blast_results_path = './data/blast/blast_results.xml'
# blast_dict = readBlastXML(blast_results_path)

# print('---------- add BLAST results to VUS_dict ----------')
# vus_dict = processor.add_blast_results(vus_dict, blast_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_path = './data/VUS_with_blast.tsv'
# processor.write_to_tsv(vus_dict, header, vus_blast_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS_with_blast.tsv'
# vus_dict = processor.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make_tsv_for_CADD ----------')
# caddfile_path = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/vus_cadd.vcf'
# processor.make_tsv_for_CADD(vus_dict, caddfile_path)
# print('Done')

# print('---------- get CADD scores from its website ----------')
# cadd_results = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/GRCh38-v1.4.tsv.gz'
# web_handler.get_CADD_scores(caddfile_path, cadd_results)

# print('---------- read CADD results ----------')
# cadd_results = './data/CADD/GRCh38-v1.4.tsv.gz'
# cadd_dict = processor.read_CADD_results(cadd_results)

# print('---------- add CADD results ----------')
# vus_dict, unfound_cadd_dict = processor.add_cadd_results(vus_dict, cadd_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_cadd_path = './data/VUS_with_blast_cadd.tsv'
# processor.write_to_tsv(vus_dict, header, vus_blast_cadd_path)

print('********** read TSV **********')
tsv_path = './data/VUS_with_blast_cadd.tsv'
vus_dict = processor.read_tsv_to_dict(tsv_path)
print(f'The number of VUS: {len(vus_dict)}')
print(vus_dict[0])

# print('---------- make_tsv_for_VEP ----------')
# vep_input_path = './data/vep/vep_vus_input.tsv'
# processor.make_tsv_for_vep(vus_dict, vep_input_path)
# print('Done')

## print('---------- run vep ----------')
## vep_output_path = './data/vep/vep_vus_results.tsv'
## vep_input_path = './data/vep/vep_vus_input_short.tsv'
## processor.vep_locally(vep_input_path, vep_output_path)

# print('---------- make ordered_dict for VEP ----------')
# vus_ordered_dict = processor.make_ordered_vus_dict_for_vep(vus_dict)

# print('---------- make_tsv_ordered_for_VEP ----------')
# vep_input_path = './data/vep/vep_vus_input_ordered.tsv'
# processor.make_tsv_ordered_for_vep(vus_ordered_dict, vep_input_path)
# print('Done')

print('---------- read VEP output ----------')
vep_output_path = './data/vep/vep_vus_output.txt'
# vep_output_path = './data/vep/vep_vus_output_short.txt'
vep_dict = processor.read_vep_output(vep_output_path)

print('---------- add VEP output ----------')
vus_dict = processor.add_vep_output(vus_dict, vep_dict)
header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
vus_blast_cadd_path = './data/VUS_with_blast_cadd_vep.tsv'
processor.write_to_tsv(vus_dict, header, vus_blast_cadd_path)
