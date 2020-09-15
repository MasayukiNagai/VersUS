#!/usr/bin/env python
# coding: utf-8

# In[35]:


import pandas as pd
import numpy as np
import csv
import os
from lxml import etree
import xml.etree.ElementTree as ET
from Bio import Entrez
from Bio import SeqIO
import time
import datetime

Entrez.email = "mnaffairs.intl@gmail.com"


# In[36]:


aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


# In[37]:


def readHumanGenesEC(path):
    genes_dict = {}
    dup = set()
    with open(path, 'r') as f:
        f.readline()
        reader = csv.reader(f)
        for row in reader:
            key = row[0]
            if key in genes_dict.keys():
                dup.add(key)
            else:
                genes_dict[key] = row[1]
    print(f'{len(dup)} genes are duplicated: ', dup)
    return genes_dict


# In[38]:


genes_dict = readHumanGenesEC('../data/HumanEnzWithEC.csv')
print(f'{len(genes_dict)} genes are imported')


# In[52]:


class variationHandler(object):
    def __init__(self, genes_dict):
        self.dictlist = {'gene_ID':[], 'gene_name':[], 'clinical_significance': [], 'EC_number': [], 'missense_variation': [], 'NP_accession': [], 'ClinVar_accession':[], 'Chr': [], 'start':[], 'stop':[], 'referenceAllele':[], 'alternateAllele':[]}
        self.genes_dict = genes_dict
        self.unnecessary_types = ('inversion', 'copy number gain', 'tandem duplication', 'microsatellite', 'copy number loss', 'distinct chromosomes', 'fusion', 'complex', 'duplication', 'translocation')
        self.is_type = False
        self.gene_ID = ""
        self.gene_name = ""
        self.clinvar_acc = ""
        self.np_acc = ""
        self.change = ""
        self.chr = ""
        self.start_num = ""
        self.stop_num = ""
        self.referenceAllele = ""
        self.alternateAllele = ""
        self.is_GeneList = False
        self.ct_gene = 0
        self.check_grch = False
        self.ct_np = 0
        self.ct_mc = 0  # counter for Molecular Consequence tag
        self.is_haplotype = False  # check if a variation is haplotype or not
        self.is_genotype = False  # check if a variation is genotype or not
        self.is_missense = False
        self.is_conflicting = False
        self.is_not_provided = False
        self.is_interpretations = False
        self.is_interpretation = False
        self.is_description = False
        self.is_desc_hist = False
        self.intpn = []
        self.ct = 0
        self.ct_missense = 0
        self.ct_uncertain = 0
        self.ct_conflicting = 0
        self.ct_not_provided = 0
        
    def start(self, tag, attrs):
        if (tag == 'VariationArchive') and (attrs.get('VariationType').lower() not in self.unnecessary_types):
            self.is_type = True
            if attrs.get('VariationType').lower() == 'haplotype':
                self.is_haplotype = True
            if attrs.get('VariationType').lower() == 'compoundheterozygote':
                self.is_genotype = True
            self.clinvar_acc = attrs.get('Accession')
        if self.is_type:
            if tag == 'GeneList':
                self.is_GeneList = True
            elif tag == 'Gene' and self.ct_gene == 0:
                self.gene_ID = attrs.get('Symbol')
                self.gene_name = attrs.get('FullName')
                self.ct_gene += 1
            elif tag == 'SequenceLocation' and self.is_GeneList == False and self.check_grch == False:
                if attrs.get('Assembly') == 'GRCh38':
                    self.chr = attrs.get('Chr')
                    self.start_num = attrs.get('start')
                    self.stop_num = attrs.get('stop')
                    self.referenceAllele = attrs.get('referenceAlleleVCF')
                    self.alternateAllele = attrs.get('alternateAlleleVCF')
                    self.check_grch = True
            elif tag == 'ProteinExpression' and self.ct_np == 0:
                self.np_acc = attrs.get('sequenceAccessionVersion')
                self.change = attrs.get('change') 
                if self.np_acc and self.np_acc.startswith('NP'):            
                    self.ct_np += 1
            elif tag == 'MolecularConsequence' and self.ct_mc == 0:
                if attrs.get('Type') and 'missense' in attrs.get('Type').lower():
                    self.is_missense = True
                    self.ct_missense += 1
                self.ct_mc += 1
            elif tag == 'Interpretations':
                self.is_interpretations = True
            elif tag == 'Interpretation':
                self.is_interpretation = True
            elif tag == 'Description':
                self.is_description = True
            elif tag == 'DescriptionHistory':
                self.is_desc_hist = True
            
    def end(self, tag):
        if (tag == 'VariationArchive' and self.is_type) or ((self.is_haplotype or self.is_genotype) and tag == 'SimpleAllele'):
            if len(self.intpn) == 1:
                clinical_significance = self.intpn[0].lower()
                if "uncertain" in clinical_significance:
                    self.is_uncertain = True
                    self.ct_uncertain += 1
                elif "conflicting" in clinical_significance:
                    self.is_conflicting = True
                    self.ct_conflicting += 1
                elif "not provided" in clinical_significance:
                    self.is_not_provided = True
                    self.ct_not_provided += 1
#             elif len(self.intpn) > 1:
#                 print(f'Interpretaion: {self.intpn}, Accession: {self.clinvar_acc}, count: {self.ct}')
            
            if (self.gene_ID in self.genes_dict.keys()) and self.is_missense and (self.is_uncertain or self.is_conflicting or self.is_not_provided):
                try:
                    self.change = self.change.split('p.')[1]
                    before = aatranlation.get(self.change[0:3])
                    after = aatranlation.get(self.change[len(self.change) - 3:len(self.change)])
                except: 
                    before = None
                    after = None
                if before and after:  # check if both have a value in aa dict
                    num = self.change[3:len(self.change) - 3]
                    abbreviated_change = before + num + after
                    fasta = np.nan
                    self.dictlist['gene_ID'].append(self.gene_ID)
                    self.dictlist['gene_name'].append(self.gene_name)
                    self.dictlist['clinical_significance'].append(clinical_significance)
                    self.dictlist['EC_number'].append(self.genes_dict.get(self.gene_ID))
                    self.dictlist['missense_variation'].append(abbreviated_change)
                    self.dictlist['NP_accession'].append(self.np_acc)
                    self.dictlist['ClinVar_accession'].append(self.clinvar_acc)
                    self.dictlist['Chr'].append(self.chr)
                    self.dictlist['start'].append(self.start_num)
                    self.dictlist['stop'].append(self.stop_num)
                    self.dictlist['referenceAllele'].append(self.referenceAllele)
                    self.dictlist['alternateAllele'].append(self.alternateAllele)
            self.ct_gene = 0             
            self.check_grch = False
            self.is_missense = False
            self.is_uncertain = False
            self.is_conflicting = False
            self.is_not_provided = False
            self.ct_np = 0
            self.ct_mc = 0
            self.intpn = []
        if self.is_type:
            if tag == 'GeneList':
                self.is_GeneList = False
            elif tag == 'Interpretations':
                self.is_interpretations = False
            elif tag == 'Interpretation':
                self.is_interpretaion = False
            elif tag == 'Description':
                self.is_description = False
            elif tag == 'DescriptionHistory':
                self.is_desc_hist = False
        if tag == 'VariationArchive':
            self.is_type = False
            self.is_haplotype = False
            self.is_genotype = False
            self.ct +=1
            if self.ct % 10000 == 0:
                print(f'counter: {self.ct}')
                
    def data(self, data):
        if self.is_interpretations and self.is_interpretation and self.is_description and (not self.is_desc_hist):
            self.intpn.append(data)
            
    def close(self):
        print(f"Variations: {self.ct}")
        print(f"Uncertain Significance: {self.ct_uncertain}")
        print(f"Conflicting Report: {self.ct_conflicting}")
        print(f"Not Provided: {self.ct_not_provided}")
        print(f"Missense: {self.ct_missense}")
        print(f"Mutations in the list: {len(self.dictlist['gene_ID'])}")
        print('debug: the file is closed')
        return self.dictlist


# In[53]:


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXML(input_path, output_path, gene_set):
    print('debug: start parcing')
    parser = etree.XMLParser(target=variationHandler(gene_set))
    data = etree.parse(input_path, parser)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index = False, header = True)
    return df


# In[55]:


xmlfile = '../data/ClinVarVariationRelease_00-latest_weekly.xml'
out_path = '../data/MM_enzyme.csv'
df_0 = readClinVarVariationsXML(xmlfile, out_path, genes_dict)


# In[54]:


# this column is a test for a certain variation I want to look over
subnode = '../data/subnode.xml'
temp_path = '../data/temp.csv'
df_temp = readClinVarVariationsXML(subnode, temp_path, genes_dict)
df_temp.head()


# In[7]:


class variationHandlerSpecific(object):
    def __init__(self, accession):
        self.is_accession = False
        self.accession = accession
        self.ct = 0
        print(self.accession)
        
    def start(self, tag, attrs):
        global WFILE
        if (tag == 'VariationArchive') and (attrs.get('Accession') == self.accession):
            self.is_accession = True
            print('The specific variation is found: ' + str(self.ct))
        if self.is_accession:
            if len(attrs.keys()) == 0:
                WFILE.write('<' + tag)
            else:   
                for i, t in enumerate(attrs.keys()):
                    if i == 0:
                        WFILE.write('<' + tag + ' ')
                    elif i != len(attrs.keys()) - 1:
                        WFILE.write(t + '="' + attrs.get(t) + '"' + ' ')
                    else:
                        WFILE.write(t + '="' + attrs.get(t) + '"')
            WFILE.write('>')
            
    def end(self, tag):
        global WFILE
        if self.is_accession:
            WFILE.write('</' + tag + '>')
        if tag == 'VariationArchive':
            if self.ct % 10000 == 0:
                print(self.ct)
            self.ct += 1
        if self.is_accession and tag == 'VariationArchive':
            self.is_accession = False
            WFILE.close()
            print('The subnode file is completed')
            
    def data(self, data):
        global WFILE
        if data is not None:
            if self.is_accession and data != "":
                WFILE.write(data)
            
    def close(self):
        print('The xml file is closed')


# In[8]:


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXMLSpecific(input_path, accession):
    print('Start parcing')
    parser = etree.XMLParser(target=variationHandlerSpecific(accession))
    etree.parse(input_path, parser)


# In[11]:


# this column is to get a certain part of the xml file
WFILE = open('../data/subnode.xml', 'w')
xmlfile = '../data/ClinVarVariationRelease_00-latest_weekly.xml'
accession = 'VCV000065613'
readClinVarVariationsXMLSpecific(xmlfile, accession)


# In[57]:


# makes dictionary of fasta sequences and np number 
# returns the dictionary
def makeDictOfFasta(dictpath):
    fasta_dict = {}
    for root, d_names, file_names in os.walk(dictpath):
        for filename in file_names:
            fname = os.path.join(root, filename)
            with open(fname, 'r') as f:
                print('opened a fasta file')
                np_num = ''
                sequence = ''
                for line in f:
                    if line[0] == '>':
                        if sequence != '':
                            fasta_dict[np_num] = sequence
                            np_num = ''
                            sequence = ''                    
                        i = 1
                        while line[i] != ' ':
                            np_num += line[i]
                            i += 1
                    else:
                        line = line.strip('\n')
                        sequence += line
    print(f'length of fasta dictionary: {len(fasta_dict)}')
    return fasta_dict


# In[58]:


fasta_dict = makeDictOfFasta('../fasta_sequences/')


# In[8]:


print(fasta_dict.get('NP_037514.2'))


# In[59]:


# crops fasta sequence
# returns the cropped sequnece with a specified range
def cropFASTA(sequence, location, reference, seqRange):
    if location - 1 < len(sequence) and sequence[location - 1] == reference:
        if location - 1 - seqRange <= 0:
            proteinSeq = sequence[0:2 * seqRange + 1]
        elif location - 1 + seqRange > len(sequence) - 1:
            proteinSeq = sequence[0 if len(sequence) - 1 - 2 * seqRange <= 0 else len(sequence) - 1 - 2 * seqRange:len(sequence)]
        else:
            proteinSeq = sequence[location - 1 - seqRange : location + seqRange]
        return proteinSeq
    else:
        return None 


# In[63]:


def addFASTAfromDict(fasta_dict, df):
    unfound = set()
    seq_list = []
    for index, row in df.iterrows():
        mutation = row['missense_variation']
        try:
            ref = mutation[0]
            try:
                location = int(mutation[1:len(mutation)-1])
            except:
                location = int(mutation.split('_')[0][1:])
            np_num = row['NP_accession']  # specify the column of np 
            sequence = fasta_dict[np_num]  # 
            seqRange = 12  # range of sequences to take
            seq = cropFASTA(sequence, location, ref, seqRange) if sequence else None
        except:
            seq = None
            unfound_np = row['NP_accession']
            unfound.add(unfound_np)
        seq_list.append(seq)
    df['FASTA_window'] = seq_list
    print(f'Unfound Sequences: {len(unfound)} {unfound}')
    return df


# In[61]:


df_0 = pd.read_csv('../data/MM_enzyme.csv')
df_0.head()


# In[64]:


df_1 = addFASTAfromDict(fasta_dict, df_0)
df_1.to_csv('../data/MM_enzyme.csv', index = False, header = True)
df_1.head()


# In[65]:


def getFASTA(ls_np):
    fasta_dict = {}
    for np_num in ls_np:
        handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
        seq_record = SeqIO.read(handle, 'fasta')
        sequence = seq_record.seq
        fasta_dict[np_num] = sequence
    return fasta_dict


# In[70]:


unfound_nps = ('NP_000114.2', 'NP_001017973.1', 'NP_000274.2','NP_000299.2', 'NP_001138503.1', 'NP_004251.3', 'NP_001276896.1', 'NP_001180240.1', 'NP_000195.2', 'NP_001017975.4', 'NP_001202.4', 'NP_945347.2', 'NP_001273003.1', 'NP_000439.1', 'NP_000450.2', 'NP_065829.3', 'NP_004597.2', 'NP_001191354.1', 'NP_001041636.1', 'NP_000342.2', 'NP_001177317.1')
add_fasta = getFASTA(unfound_nps)
fasta_dict.update(add_fasta)
      


# In[71]:


df_1 = addFASTAfromDict(fasta_dict, df_0)
df_1.to_csv('../data/MM_enzyme.csv', index = False, header = True)
df_1.head()


# In[84]:


df_1[['NP_accession', 'gene_ID', 'gene_name', 'FASTA_window']].head()


# In[85]:


# make FASTA format text file from dataframe for blast search
def makeFASTAfile(df, output_path):
    subset = df[['NP_accession', 'gene_ID', 'gene_name', 'FASTA_window']]
    tuples = [tuple(x) for x in subset.values]
    with open(output_path, 'w') as f:
        for tup in tuples:
            line = '>' + '\t'.join(tup[0:3]) + '\n'
            fasta = str(tup[3]) + '\n'
            f.write(line)
            f.write(fasta) 


# In[86]:


makeFASTAfile(df_1, '../data/fasta.txt')


# In[87]:


def blastLocal(fasta_path, out_path, evalue=10.0, window_size=3):
    cmd1 = '../ncbi/blast/bin/'
    cmd2 = 'blastp' + ' '         + '-query ' + fasta_path + ' '         + '-db ' + cmd1 + 'pdbaa' + ' '         + '-evalue ' + str(evalue) + ' '         + '-outfmt ' + '5' + ' '         + '-out ' + out_path
    cmd = cmd1 + cmd2
    b_cmd = os.system(cmd)
    print(cmd + ' : ran with exit code %d' %b_cmd)


# In[88]:


start = datetime.datetime.now()

fasta_path = '../data/fasta.txt'  # path for the fasta file
out_path = '../data/blast_result.xml'  # path for the output file
evalue = 10.0
size = 3
blastLocal(fasta_path, out_path, evalue, size)

end = datetime.datetime.now()
time = end - start
c = divmod(time.days * 86400 + time.seconds, 60)
print(c)


# In[93]:


class BlastHandler(object):
    def __init__(self):
        self.dictlist = {'pdb_ID': [], 'BLAST_evalue': [], 'hit_from': [], 'hit_to': []}
        self.is_hit_id = False
        self.is_evalue = False
        self.is_hit_from = False
        self.is_hit_to = False
        self.ct_iter = 0
        self.ct = 0
        
    def start(self, tag, attrs):
        if tag == 'Hit':  # there are a few <Hsp> tags within a <Hit> tag
            self.ct_iter += 1
        elif tag == 'Hit_id':
            self.is_hit_id = True
        elif tag == 'Hsp_evalue':
            self.is_evalue = True
        elif tag == 'Hsp_hit-from':
            self.is_hit_from = True
        elif tag == 'Hsp_hit-to':
            self.is_hit_to = True
            
    def end(self, tag):
        if tag == 'Iteration':
            if self.ct_iter == 0:  # when there is no hit, add None to each list
                self.dictlist['pdb_ID'].append(None)
                self.dictlist['BLAST_evalue'].append(None)
                self.dictlist['hit_from'].append(None)
                self.dictlist['hit_to'].append(None)
            self.ct_iter = 0
            self.ct += 1
            if self.ct % 10000 == 0:
                print(f'coutner: {self.ct}')
        elif tag == 'Hit_id':
            self.is_hit_id = False
        elif tag == 'Hsp_evalue':
            self.is_evalue = False
        elif tag == 'Hsp_hit-from':
            self.is_hit_from = False
        elif tag == 'Hsp_hit-to':
            self.is_hit_to = False
        elif tag == 'Hsp':
            self.ct_iter += 1
        
    def data(self, data):
        if self.ct_iter == 1:
            if self.is_hit_id:
                try:
                    pdb = data.split('pdb|')[1]
                except: 
                    print(f'Cannot split pdb {data}')
                    pdb = None
                self.dictlist['pdb_ID'].append(pdb)
            elif self.is_evalue:
                self.dictlist['BLAST_evalue'].append(data)
            elif self.is_hit_from:
                self.dictlist['hit_from'].append(data)
            elif self.is_hit_to:
                self.dictlist['hit_to'].append(data)
    
    def close(self):
        print('Has completed parsing the blast_result xml')
        print(f'Total Count: {self.ct}')
        print(f'Length pdb: {len(self.dictlist["pdb_ID"])}, evalue: {len(self.dictlist["BLAST_evalue"])}, hit_from : {len(self.dictlist["hit_from"])}, hit_to : {len(self.dictlist["hit_to"])}')
        return self.dictlist


# In[94]:


def readBlastXML(input_path):
    print('Start parcing')
    parser = etree.XMLParser(target=BlastHandler())
    data = etree.parse(input_path, parser)
    return data    


# In[95]:


def combineBlastResults(df, dictlist):
    df['pdb_ID'] = dictlist['pdb_ID']
    df['BLAST_evalue'] = dictlist['BLAST_evalue']
    df['hit_from'] = dictlist['hit_from']
    df['hit_to'] = dictlist['hit_to']
    return df


# In[96]:


blast_dl = readBlastXML('../data/blast_result.xml')
blast_df = pd.DataFrame(blast_dl)
print(len(blast_dl['pdb_ID']))
blast_df.head()


# In[97]:


df_1 = pd.read_csv('../data/MM_enzyme.csv')
df_2 = combineBlastResults(df_1, blast_dl)
df_2.to_csv('../data/MM_enzyme_blast.csv', index = False, header = True)
df_2.head()


# In[98]:


# df_2 = pd.read_csv('../data/MM_enzyme_blast.csv')
# df_2[df_2.duplicated(['ClinVar_accession'])]


# In[126]:


columns = ('gene_ID', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'SIFT_prediction', 'CADD_score', 'Chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
columns2 = ('gene_name', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'Chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
df_2.loc[:, columns2].head()   


# In[127]:


def orderColumns(df, col_names):
    df_cols = set(df.colums)
    if set(col_names).issubset(df_cols):
            df_ordered = df.loc[:, col_names]
            return df_ordered
    else:
        print('column names given include missing labels')


# In[101]:


