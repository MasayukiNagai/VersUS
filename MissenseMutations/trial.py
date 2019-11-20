import pandas as pd
import numpy as np
import csv
import xml.etree.ElementTree as ET
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW

Entrez.email = "mnaffairs.intl@gmail.com"

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


def getHumanEnzymeIDs():
    handleHumanEnzymes = Entrez.esearch(db="protein", retmax=20000,
                                        term="Homo sapiens[Organism] AND RefSeq[Filter] AND (1*[EC/RN Number] OR 2*[EC/RN Number] OR 3*[EC/RN Number] OR 4*[EC/RN Number] OR 5*[EC/RN Number] OR 6*[EC/RN Number])")
    readHumanEnzymes = Entrez.read(handleHumanEnzymes)
    return readHumanEnzymes['IdList']


def getHumanGenes(idLists, path):
    human_enzymes = set()
    protein_info = Entrez.efetch(db='protein', id=idLists, retmode='xml', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
    tree = ET.parse(protein_info)
    root = tree.getroot()
    count = 0  # debug
    for gene in root.findall('./GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier/[GBQualifier_name="gene"]/GBQualifier_value'):
        human_enzymes.add(gene.text)
        count += 1  # debug
    print(count)
    with open('../data/HumanEnzymes.txt', 'w') as filehandle:
        filehandle.seek(0)  # sets the file's current position at 0
        filehandle.truncate()  # sets the logical size of a file to 0, which means to delte its content
        filehandle.writelines('%s\n' % human_enzyme for human_enzyme in human_enzymes)
    return human_enzymes


def readHumanGenes(path):
    human_genes = []
    with open(path, 'r') as filehandle:
        human_genes = filehandle.read().splitlines()
    return human_genes


def getVUSIDs():
    handleVUS = Entrez.esearch(db='clinvar', retmax=200000,
                                    term='( ( ("clinsig has conflicts"[Properties]) OR ("clinsig vus"[Properties]) ) AND ("missense variant"[molecular consequence] OR "SO 0001583"[molecular consequence]))')
    recordVUS = Entrez.read(handleVUS)
    return recordVUS['IdList']


def getVUS_Clinvar(idLists, genes, path):
    data = {'Gene': [], 'VUS_protein': [], 'NP': [], 'PDB': []}
    i = 0
    count = 0  # debug
    geneSet = set()
    while i < len(idLists):
        missense_info = Entrez.efetch(db='clinvar', id=idLists[i:i+10000], retmax=10000, rettype='vcv', is_variationid="true", from_esearch="true",
                                  api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
        tree = ET.parse(missense_info)
        root = tree.getroot()
        for variation in root.findall('./VariationArchive'):
            gene = variation.find('./InterpretedRecord/SimpleAllele/GeneList/Gene')
            name = gene.attrib['Symbol']
            count += 1
            if name in genes:
                geneSet.add(name)
                for mutation in variation.findall('.//ProteinExpression'):
                    np_num = mutation.attrib['sequenceAccessionVersion']
                    change = mutation.attrib['change'].split('p.')[1]
                    before = aatranlation.get(change[0:3])
                    after = aatranlation.get(change[len(change) - 3:len(change)])
                    if before and after:  # check if both have a value in aa dict
                        num = change[3:len(change) - 3]
                        abbreviated_change = before + num + after
                        data['Gene'].append(name)
                        data['VUS_protein'].append(abbreviated_change)
                        data['NP'].append(np_num)
                        data['PDB'].append(np.nan)
        if i + 10000 < len(idLists):
            i += 10000
        elif i == len(idLists) - 1:
            break
        else:
            i += len(idLists) % 10000 - 1
        print(count)  # debug
    print(len(geneSet))                
    df = pd.DataFrame(data)
    df.to_csv(path, index = False, header = True)
    return df


def readVUScsv(path):
    data = []
    with open(path) as filehandle:
        reader = csv.reader(filehandle)
        data = list(reader)
    df = pd.DataFrame(data)
    return df


def getFASTA(np_num, location, beforeMutation, numOfSequence = 10):
    handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text')
    seq_record = SeqIO.read(handle, 'fasta')
    sequence = seq_record.seq
    if location - 1 < len(sequence) and sequence[location - 1] == beforeMutation:
        proteinSeq = sequence[0 if location - 1 - numOfSequence <= 0 else location - 1 - numOfSequence : location + numOfSequence]
        return proteinSeq
    else:
        print("error")  # debug
        return False

id = 558835
missense_info = Entrez.efetch(db='clinvar', id=id, retmax=10000, rettype='vcv', is_variationid="true", from_esearch="true",
                              api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
readHumanGenes = missense_info.read()
print(readHumanGenes)

# sequence = getFASTA('NP_005557.1', 190, 'L', 10)
# print(sequence)

# handle = NCBIWWW.qblast("blastp", 'pdb', proteinSeq, expect = 10.0)
# Hit_id