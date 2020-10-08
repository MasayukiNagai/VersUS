# test script
from main import *
from XMLHandler import *
from Processor import * 


processor = Processor()

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
# caddfile_path = './data/CADD/vus_cadd.vcf'
# processor.make_tsv_for_CADD(vus_dict, caddfile_path)
# print('Done')

# print('---------- read CADD results ----------')
# cadd_results = './data/CADD/GRCh38-v1.4.tsv'
# cadd_dict = processor.read_CADD_results(cadd_results)

# print('---------- add CADD results ----------')
# vus_dict, unfound_cadd_dict = processor.add_cadd_results(vus_dict, cadd_dict)
# print('Unfound cadd: ', unfound_cadd_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_cadd_path = './data/VUS_with_blast_cadd.tsv'
# processor.write_to_tsv(vus_dict, header, vus_blast_cadd_path)