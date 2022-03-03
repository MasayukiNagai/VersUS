# test script
import os
import logging
from Handlers.AlphaFoldHandler import AlphaFoldHandler
from Handlers.SeqHandler import *
from Handlers.ClinVarHandler import *
from Handlers.BLASTHandler import *
from Handlers.CADDHandler import *
from Handlers.VEPHandler import *
from Handlers.FileHandler import *
from VersUS import *

gene_file = './data/gene/HumanEnzWithEC.csv'
seq_dir = './data/fasta_sequences'
blast_path = './ncbi/blast/bin'
blast_input = './data/blast/vus_blast.fasta'
blast_output = './data/blast/blast_results.xml'
cadd_input = '/Users/moon/DePauw/ITAP/VersUS/src/data/CADD/cadd_input.vcf.gz'
cadd_output = '/Users/moon/DePauw/ITAP/VersUS/src/data/CADD/CADD_scores.tar.gz'
vep_input = './data/vep/vep_vus_input.tsv'
vep_output = './data/vep/vep_vus_results.tsv'

# seqHandler = SeqHandler(gene_file, seq_dir)
# clinvarHandler = ClinVarHandler()
# blastHandler = BLASTHandler(blast_input, blast_output)
# caddHandler = CADDHandler(cadd_input, cadd_output)
# vepHandler = VEPHandler(vep_input, vep_output)
versus = VersUS()

# print('---------- read Gene File ----------')
# genes_dict = seqHandler.readHumanGenesEC()

# print('---------- VariationParser ----------')
# clinvar_file = './data/clinvar/ClinVarVariationRelease_00-latest_weekly.xml'
# # clinvar_file = './data/clinvar/clinvar_tail_5.xml'
# vus_dict = clinvarHandler.readClinVarVariationsXML(clinvar_file, genes_dict)

# print('---------- write to TSV ----------')
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele')
# outfile_path = './data/VUS.tsv'
# versus.write_to_tsv(vus_dict, header, outfile_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS.tsv'
# vus_dict = versus.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make sequence dict ----------')
# seq_dir_path = './data/fasta_sequences'
# seq_dict = seqHandler.make_seq_dict()
# print(next(iter(seq_dict.keys())))

# print('---------- add seq to vus_dict ----------')
# seq_range = 12
# vus_dict, unfound_seq = seqHandler.add_seq_to_dict(vus_dict, seq_range)
# print(f'NP numbers missing their seq: {(unfound_seq)}')

# print('---------- get unfound sequences ----------')
# seq_dict = seqHandler.fetch_seq(unfound_seq)
# vus_dict, unfound_seq = seqHandler.add_seq_to_dict(vus_dict, seq_dict)
# print(f'NP numbers missing their seq: {(unfound_seq)}')

# print('---------- write to TSV ----------')
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window')
# outfile_path = './data/VUS.tsv'
# versus.write_to_tsv(vus_dict, header, outfile_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS.tsv'
# vus_dict = versus.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make a fasta file for blast ----------')
# blastHandler.make_fasta_for_blast(vus_dict)
# print('Done')

# print('---------- run BLAST ----------')
# blastHandler.blast_locally()

# print('---------- read BLAST results ----------')
# blast_dict = blastHandler.readBlastXML()

# print('---------- add BLAST results to VUS_dict ----------')
# vus_dict = blastHandler.add_blast_results(vus_dict, blast_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_path = './data/VUS_with_blast.tsv'
# versus.write_to_tsv(vus_dict, header, vus_blast_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS_with_blast.tsv'
# vus_dict = read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make_tsv_for_CADD ----------')
# caddHandler.make_tsv_for_CADD(vus_dict)
# print('Done')

# print('---------- get CADD scores from its website ----------')
# cadd_results = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/GRCh38-v1.4.tsv.gz'
# caddHandler.get_CADD_scores()

# print('---------- read CADD results ----------')
# cadd_results = './data/CADD/GRCh38-v1.4.tsv.gz'
# cadd_dict = caddHandler.read_CADD_results()

# print('---------- add CADD results ----------')
# vus_dict = caddHandler.add_cadd_results(vus_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_cadd_path = './data/VUS_with_blast_cadd.tsv'
# versus.write_to_tsv(vus_dict, header, vus_blast_cadd_path)

# print('********** read TSV **********')
# tsv_path = './data/VUS_with_blast_cadd.tsv'
# vus_dict = versus.read_tsv_to_dict(tsv_path)
# print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

# print('---------- make_tsv_for_VEP ----------')
# vepHandler.make_tsv_for_vep(vus_dict)
# print('Done')

# print('---------- run vep ----------')
# vepHandler.vep_locally()

# print('---------- make ordered_dict for VEP ----------')
# vus_ordered_dict = vepHandler.make_ordered_vus_dict_for_vep(vus_dict)

# print('---------- make_tsv_ordered_for_VEP ----------')
# vepHandler.make_tsv_ordered_for_vep(vus_ordered_dict)
# print('Done')

# print('---------- read VEP output ----------')
# vep_dict = vepHandler.read_vep_output()

# print('---------- add VEP output ----------')
# vus_dict = vepHandler.add_vep_output(vus_dict)
# header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
# vus_blast_cadd_path = './data/VUS_with_blast_cadd_vep.tsv'
# versus.write_to_tsv(vus_dict, header, vus_blast_cadd_path)

print('********** read TSV **********')
tsv_path = '../results/vus-06132021.tsv'
vus_dict = read_tsv_to_dict(tsv_path)
print(f'The number of VUS: {len(vus_dict)}')
# print(vus_dict[0])

print('---------- AlphaFold Link Check ----------')
afhandler = AlphaFoldHandler()
vus_dict = afhandler.run(vus_dict)
header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to', 'alphafold')
vus_output = './data/VersUS_alphafold.tsv'
versus.write_to_tsv(vus_dict, header, vus_output)
