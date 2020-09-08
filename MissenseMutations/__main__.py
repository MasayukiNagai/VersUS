import pandas as pd
import numpy as np
import csv
import os
import xml.etree.ElementTree as ET
from lxml import etree
from Bio import Entrez
from Bio import SeqIO
# from Bio.Blast import NCBIWWW
import time
import datetime

Entrez.email = "mnaffairs.intl@gmail.com"

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


# read csv file of gene names of human enzymes
# return dict of gene names
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
    print(f'{len(dup)} genes are duplicated: ', dup)  # prints, if any, genes duplicated in the file
    return genes_dict


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


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXML(input_path, output_path, gene_set):
    print('debug: start parcing')
    parser = etree.XMLParser(target=variationHandler(gene_set))
    data = etree.parse(input_path, parser)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index = False, header = True)
    return df


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


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXMLSpecific(input_path, accession):
    print('Start parcing')
    parser = etree.XMLParser(target=variationHandlerSpecific(accession))
    etree.parse(input_path, parser)

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


# fetches fasta sequences for varinants whose sequences aren't in the imported files
# returns dict of fasta
def getFASTA(ls_np):
    fasta_dict = {}
    for np_num in ls_np:
        handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
        seq_record = SeqIO.read(handle, 'fasta')
        sequence = seq_record.seq
        fasta_dict[np_num] = sequence
    return fasta_dict


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


# make FASTA format text file from dataframe for blast search
def makeFastaFileForBlast(df, output_path):
    subset = df[['NP_accession', 'gene_ID', 'gene_name', 'FASTA_window']]
    tuples = [tuple(x) for x in subset.values]
    with open(output_path, 'w') as f:
        for tup in tuples:
            line = '>' + '\t'.join(tup[0:3]) + '\n'
            fasta = str(tup[3]) + '\n'
            f.write(line)
            f.write(fasta) 


def blastLocally(fasta_path, out_path, evalue=10.0, window_size=3):
    cmd1 = '../ncbi/blast/bin/'
    cmd2 = 'blastp' + ' '         + '-query ' + fasta_path + ' '         + '-db ' + cmd1 + 'pdbaa' + ' '         + '-evalue ' + str(evalue) + ' '         + '-outfmt ' + '5' + ' '         + '-out ' + out_path
    cmd = cmd1 + cmd2
    b_cmd = os.system(cmd)
    print(cmd + ' : ran with exit code %d' %b_cmd)


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


def readBlastXML(input_path):
    print('Start parcing')
    parser = etree.XMLParser(target=BlastHandler())
    data = etree.parse(input_path, parser)
    return data    


def addBlastResults(df, dictlist):
    df['pdb_ID'] = dictlist['pdb_ID']
    df['BLAST_evalue'] = dictlist['BLAST_evalue']
    df['hit_from'] = dictlist['hit_from']
    df['hit_to'] = dictlist['hit_to']
    return df


def orderColumns(df, col_names):
    df_cols = set(df.colums)
    if set(col_names).issubset(df_cols):
            df_ordered = df.loc[:, col_names]
            return df_ordered
    else:
        print('column names given include missing labels')


# read csv file
# return dataframe
def readCSV(path):
    return pd.read_csv(path)






def removeSameVariations(path):
    data = []
    with open(path) as filehandle:
        reader = csv.reader(filehandle)
        data = list(reader)
    length = len(data)
    ct = 0
    i = 1
    while ct < length - 1:
        geneCurr = data[i][0]
        geneNext = data[i+1][0]
        if geneCurr == geneNext:
            mutation = data[i][1]
            before = mutation[0]
            afterCurr = mutation[len(mutation)-1]
            location = int(mutation[1:len(mutation)-1])
            sequenceCurr = getFASTA(data[i][2], location, before, 10)
            mutation = data[i+1][1]
            before = mutation[0]
            afterNext = mutation[len(mutation)-1]
            location = int(mutation[1:len(mutation)-1])
            sequenceNext = getFASTA(data[i+1][2], location, before, 10)
            if sequenceCurr == sequenceNext and afterCurr == afterNext:
                data.pop(i+1)
                print(geneCurr)
                print("popped")
            else:
                i += 1
        else:
            i += 1
        ct += 1
        print(ct)
    df = pd.DataFrame(data)
    df.to_csv(path, index = False, header= False)
    return df


def removeSameVariationsByName(inpath, outpath):
    with open(outpath, 'w') as outFile, open(inpath, 'r') as inputFile:
        writer = csv.writer(outFile, delimiter=',')
        reader = csv.reader(inputFile, delimiter=',')
        writer.writerow(next(reader))
        rowCurr = next(reader)
        print(rowCurr)
        for rowNext in reader:
            geneCurr = rowCurr[0]
            geneNext = rowNext[0]
            np = rowCurr[2]
            print(geneCurr + " " + geneNext)
            if geneCurr == geneNext:
                mutationCurr = rowCurr[1]
                mutationNext = rowNext[1]
                if mutationCurr[0:len(mutationCurr)-2] != mutationNext[0:len(mutationNext)-2]:
                    if np[0:3] == 'NP_':
                        writer.writerow(rowCurr)
                        print("debug: added")
                        rowCurr = rowNext
            else:
                if np[0:3] == 'NP_':
                    writer.writerow(rowCurr)
                rowCurr = rowNext


def removeSameVariationsBySequence(inpath, outpath):
    with open(outpath, 'w') as outFile, open(inpath, 'r') as inputFile:
        writer = csv.writer(outFile, delimiter=',')
        reader = csv.reader(inputFile, delimiter=',')
        writer.writerow(next(reader))
        rowCurr = next(reader)
        sequenceCurr = ''
        print(rowCurr)
        for rowNext in reader:
            geneCurr = rowCurr[0]
            geneNext = rowNext[0]
            print(geneCurr + " " + geneNext)
            if geneCurr == geneNext:
                mutationCurr = rowCurr[1]
                beforeCurr = mutationCurr[0]
                afterCurr = mutationCurr[len(mutationCurr)-1]
                locationCurr = int(mutationCurr[1:len(mutationCurr)-1])
                npCurr = rowCurr[2]

                mutationNext = rowNext[1]
                beforeNext = mutationNext[0]
                afterNext = mutationNext[len(mutationNext)-1]
                locationNext = int(mutationNext[1:len(mutationNext)-1])
                npNext = rowNext[2]
                if beforeCurr == beforeNext:
                    sequenceNext = getFASTA(npNext, locationNext, beforeNext, 10)
                    if sequenceCurr != sequenceNext:
                        writer.writerow(rowCurr)
                        print("debug: added")
                        rowCurr = rowNext
                        sequenceCurr = sequenceNext
                else:
                    print("debug: added")
                    writer.writerow(rowCurr)
                    rowCurr = rowNext
            else:
                rowCurr = rowNext


    # data = list(reader)
    # length = len(data)
    # ct = 0
    # i = 1
    # while ct < length - 1:
    #     geneCurr = data[i][0]
    #     geneNext = data[i+1][0]
    #     if geneCurr == geneNext:
    #         mutation = data[i][1]
    #         before = mutation[0]
    #         afterCurr = mutation[len(mutation)-1]
    #         location = int(mutation[1:len(mutation)-1])
    #         sequenceCurr = getFASTA(data[i][2], location, before, 10)
    #         mutation = data[i+1][1]
    #         before = mutation[0]
    #         afterNext = mutation[len(mutation)-1]
    #         location = int(mutation[1:len(mutation)-1])
    #         sequenceNext = getFASTA(data[i+1][2], location, before, 10)
    #         if sequenceCurr == sequenceNext and afterCurr == afterNext:
    #             data.pop(i+1)
    #             print(geneCurr)
    #             print("popped")
    #         else:
    #             i += 1
    #     else:
    #         i += 1
    #     ct += 1
    #     print(ct)
    # df = pd.DataFrame(data)
    # df.to_csv(path, index = False, header= False)
    # return df                    


# get every enzyme
# enzyme_ids = getHumanEnzymeIDs()
# print(str(len(enzyme_ids)) + " enzymes are found".format(len(enzyme_ids)))


# get human genes
# human_genes = getHumanGenes(enzyme_ids, '../data/HumanEnzymes.txt')
# print(human_genes)
# print(len(human_genes))


# read human genes text file
# human_genes = readHumanGenes('../data/UniProtHumanEnzymeGenes.txt')
# print(human_genes)
# print(str(len(human_genes)) + " genes of human enzymes are imported")


# get every missense genes
# VUS_ids = getVUSIDs()
# print(str(len(VUS_ids)) + " variants are found")


# get csv file which filters VUS_ids out with human_genes  
# df_VUS = fetchVUS_ClinVar(VUS_ids, human_genes, '../data/MM_enzyme.csv')
# print(df_VUS)

# addPDB('../data/MM_enzyme.csv')

# df = removeSameVariations('../data/MM_enzyme.csv')
# print(df)

# removeSameVariationsByName('../data/MM_enzyme.csv', '../data/MM_enzyme_filter_1.csv')

# removeSameVariationsBySequence('../data/MM_enzyme_filter_1.csv', '../data/MM_enzyme_filter_2.csv')
# start = datetime.datetime.now()
# sequence = getFASTA('NP_005557.1', 190, 'L', 10)
# sequence = 'KFGELVAEEARRKGELRYMHS'
# pdb = getPDB(sequence, 10.0)
# print(pdb)
# end = datetime.datetime.now()
# time = end - start
# c = divmod(time.days * 86400 + time.seconds, 60)
# print(c)