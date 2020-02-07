import csv 
import pandas as pd
import numpy as np
from lxml import etree

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

# read csv file of gene names of human enzymes
# return list of the enzymes
def readHumanGenes(path):
    human_genes = []
    with open(path, 'r') as filehandle:
        human_genes = filehandle.read().splitlines()
    return human_genes


class variationHandler(object):
    def __init__(self, enzyme_genes):
        self.dictlist = {'interpretation': [], 'gene':[], 'accession':[], 'mutation': [], 'NP': [], 'Chr': [], 'start':[], 'stop':[], 'referenceAllele':[], 'alternateAllele':[], 'FASTA':[], 'PDB': []}
        self.enzyme_genes = enzyme_genes
        self.gene = ""
        self.interpretation = ""
        self.accession = ""
        self.mutation = ""
        self.np_num = ""
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
        self.is_missense = False
        self.ct = 0

    def start(self, tag, attrs):
        if tag == 'VariationArchive':
            self.accession = attrs['Accession']
        elif tag == 'GeneList':
            self.is_GeneList = True
        elif tag == 'Gene' and self.ct_gene == 0:
            if attrs['Symbol'] in self.enzyme_genes:
                self.gene = attrs.get('Symbol')
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
            self.np_num = attrs.get('sequenceAccessionVersion')
            self.change = attrs.get('change') 
            if self.np_num.startswith('NP'):            
                self.ct_np += 1
        elif tag == 'MolecularConsequence':
            self.is_missense = True if attrs.get('Type') == 'missense variant' else False
        elif tag == 'RCVAccession':
            self.interpretation = attrs.get('Interpretation')

    def end(self, tag):
        if tag == 'VariationArchive':
            if (self.gene in self.enzyme_genes) and self.is_missense and (("Uncertain" in self.interpretation) or ("Conflicting" in self.interpretation)):
                try:
                    self.change= self.change.split('p.')[1]
                    before = aatranlation.get(self.change[0:3])
                    after = aatranlation.get(self.change[len(self.change) - 3:len(self.change)])
                except: 
                    before = None
                    after = None
                if before and after:  # check if both have a value in aa dict
                    num = self.change[3:len(self.change) - 3]
                    abbreviated_change = before + num + after
                    # fasta = getFASTA(np_num, int(num) ,before)
                    fasta = np.nan
                    self.dictlist['interpretation'].append(self.interpretation)
                    self.dictlist['gene'].append(self.gene)
                    self.dictlist['accession'].append(self.accession)
                    self.dictlist['mutation'].append(abbreviated_change)
                    self.dictlist['NP'].append(self.np_num)
                    self.dictlist['Chr'].append(self.chr)
                    self.dictlist['start'].append(self.start_num)
                    self.dictlist['stop'].append(self.stop_num)
                    self.dictlist['referenceAllele'].append(self.referenceAllele)
                    self.dictlist['alternateAllele'].append(self.alternateAllele)
                    self.dictlist['FASTA'].append(fasta)
                    self.dictlist['PDB'].append(np.nan)     
            self.check_grch = False 
            self.is_missense = False
            self.ct_np = 0
            self.ct +=1 
            if self.ct % 10000 == 0:
                print(self.ct)
        elif tag == 'GeneList':
            self.is_GeneList = False
            self.ct_gene = 0
            
    def close(self):
        print('debug: the file is closed')
        return self.dictlist


# read xml file including name, mutation, chromosome, np number
# return dataframe and write to a csv file
def readVUS_ClinVar(input_path, output_path, gene_set):
    print('debug: start parcing')
    parser = etree.XMLParser(target=variationHandler(gene_set))
    data = etree.parse(input_path, parser)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index = False, header = True)
    return df

# read human genes text file
human_genes = readHumanGenes('../data/UniProtHumanEnzymeGenes.txt')
print(str(len(human_genes)) + " genes of human enzymes are imported")

readVUS_ClinVar('../data/ClinVarVariationRelease_00-latest_weekly.xml', '../data/MM_enzyme.csv', human_genes)
