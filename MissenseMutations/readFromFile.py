import pandas as pd
import numpy as np
import csv
import xml.etree.ElementTree as ET
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import time
import datetime

Entrez.email = "mnaffairs.intl@gmail.com"

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


# fetch id numbers of human enzymes from ncbi protein database
# return list of the ids
def getHumanEnzymeIDs():
    handleHumanEnzymes = Entrez.esearch(db="protein", retmax=20000,
                                        term="Homo sapiens[Organism] AND RefSeq[Filter] AND (1*[EC/RN Number] OR 2*[EC/RN Number] OR 3*[EC/RN Number] OR 4*[EC/RN Number] OR 5*[EC/RN Number] OR 6*[EC/RN Number])")
    readHumanEnzymes = Entrez.read(handleHumanEnzymes)
    return readHumanEnzymes['IdList']


# fetch gene names from ncbi protein databse using id lists of the enzymes
# return set of the enzymes and write to a csv file
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


# read csv file of gene names of human enzymes
# return list of the enzymes
def readHumanGenes(path):
    human_genes = []
    with open(path, 'r') as filehandle:
        human_genes = filehandle.read().splitlines()
    return human_genes


# read xml file including name, mutation, chromosome, np number
# return dataframe and write to a csv file
def readVUS_ClinVar(input_path, output_path, gene_set):
    data = {'interpretation': [], 'gene':[], 'accession':[], 'mutation': [], 'NP': [], 'Chr': [], 'start':[], 'stop':[], 'referenceAllele':[], 'alternateAllele':[], 'FASTA':[], 'PDB': []}
    i = 0
    count = 0  # debug
    tree = ET.parse(input_path)
    root = tree.getroot()
    for variation in root.findall('./VariationArchive'):
        accession = variation.attrib['Accession']
        gene = variation.find('./InterpretedRecord/SimpleAllele/GeneList/Gene').attrib['Symbol']
        count += 1  # debug
        if gene in gene_set:
            for HGVS in variation.findall('.//HGVS'):
                if HGVS.attrib['Type'] == 'coding' and HGVS.find('./MolecularConsequence').attrib['Type'] == 'missense variant':
                    np_num = HGVS.find('./ProteinExpression').attrib['sequenceAccessionVersion']
                    change = HGVS.find('./ProteinExpression').attrib['change'].split('p.')[1] 
                    interpretation = variation.find('.//RCVAccession').attrib['Interpretation']
                    if ("Uncertain" in interpretation) or ("Conflicting" in interpretation):  # Uncertain and Conflicting can be full words
                        sequenceLocation = variation.find('./InterpretedRecord/SimpleAllele/Location/SequenceLocation')
                        Chr = sequenceLocation.attrib['Chr']
                        start = sequenceLocation.attrib['start']
                        stop = sequenceLocation.attrib['stop']
                        referenceAllele = sequenceLocation.attrib['referenceAlleleVCF']
                        alternateAllele = sequenceLocation.attrib['alternateAlleleVCF']
                        before = aatranlation.get(change[0:3])
                        after = aatranlation.get(change[len(change) - 3:len(change)])
                        if before and after:  # check if both have a value in aa dict
                            num = change[3:len(change) - 3]
                            abbreviated_change = before + num + after
                            # fasta = getFASTA(np_num, int(num) ,before)
                            fasta = np.nan
                            data['interpretation'].append(interpretation)
                            data['gene'].append(gene)
                            data['accession'].append(accession)
                            data['mutation'].append(abbreviated_change)
                            data['NP'].append(np_num)
                            data['Chr'].append(Chr)
                            data['start'].append(start)
                            data['stop'].append(stop)
                            data['referenceAllele'].append(referenceAllele)
                            data['alternateAllele'].append(alternateAllele)
                            data['FASTA'].append(fasta)
                            data['PDB'].append(np.nan)
                            break
        print(f'{count}, {accession}')  # debug
        if i % 10000 == 0:
            df_temp = pd.DataFrame(data)
            df_temp.to_csv(output_path, index = False, header = True)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index = False, header = True)
    return df


# read csv file of missense information
# return dataframe
def readVUScsv(path):
    data = []
    with open(path) as filehandle:
        reader = csv.reader(filehandle)
        data = list(reader)
    df = pd.DataFrame(data)
    return df


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


# ideally add fasta sequecne to csv file but not working right now
def addFASTA(path):
    data = []
    with open(path) as filehandle:
        reader = csv.reader(filehandle)
        data = list(reader)
    i = 1
    while i < len(data):
        for j in range(0, 500): 
            mutation = data[i+j][1]
            before = mutation[0]
            location = int(mutation[1:len(mutation)-1])
            sequence = getFASTA(data[i+j][2], location, before, 10)
            data[i+j][3] = sequence
            print("debug: " + str(i+j))
        print("debug: waiting " + str(i+j))
        df = pd.DataFrame(data)
        df.to_csv(path, index = False, header = False)
        i += 500
        time.sleep(5)
    # for i in range(1, len(data)):
    #     mutation = data[i][1]
    #     before = mutation[0]
    #     location = int(mutation[1:len(mutation)-1])
    #     sequence = getFASTA(data[i][2], location, before, 10)
    #     data[i][3] = sequence
    #     print("debug:" + str(i))
    df = pd.DataFrame(data)
    df.to_csv(path, index = False, header = False)
    return df


# fetch specified range of fasta sequence
# return fasta sequence if found or None
def getFASTA(np_num, location, beforeMutation, numOfSequence = 10):
    handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
    seq_record = SeqIO.read(handle, 'fasta')
    sequence = seq_record.seq
    if location - 1 < len(sequence) and sequence[location - 1] == beforeMutation:
        proteinSeq = sequence[0 if location - 1 - numOfSequence <= 0 else location - 1 - numOfSequence : location + numOfSequence]
        return proteinSeq
    else:
        return None


# blast pdb database
# return id of the pdb sequence if found or None
def getPDB(sequence, expectValue):
    pdb = None
    print("start blast")
    handle = NCBIWWW.qblast("blastp", 'pdb', sequence, expect = expectValue, hitlist_size=3)
    record = handle.read()
    root = ET.fromstring(record)
    hit_id = root.find("./BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_id")
    if hit_id is not None:
        pdb = hit_id.text.split("|")[3]
    handle.close()
    return pdb


# blast pdb database and add id of the pdb to dataframe imported from csv file
# return dataframe and write to a new csv file
def addPDB(path):
    data = []
    with open(path) as filehandle:
        reader = csv.reader(filehandle)
        data = list(reader)   
    for i in range(1,len(data)):  # try smaller range 
        mutation = data[i][1]
        before = mutation[0]
        location = int(mutation[1:len(mutation)-1])
        sequence = getFASTA(data[i][2], location, before, 10)
        print(sequence)
        data[i][4] = getPDB(sequence, 10.0)
        # data[i][3] = i
    df = pd.DataFrame(data)
    df.to_csv(path, index = False, header = False)
    return df



# get every enzyme
# enzyme_ids = getHumanEnzymeIDs()
# print(str(len(enzyme_ids)) + " enzymes are found".format(len(enzyme_ids)))


# get human genes
# human_genes = getHumanGenes(enzyme_ids, '../data/HumanEnzymes.txt')
# print(human_genes)
# print(len(human_genes))


# read human genes text file
human_genes = readHumanGenes('../data/UniProtHumanEnzymeGenes.txt')
print(str(len(human_genes)) + " genes of human enzymes are imported")


# get every missense genes
# VUS_ids = getVUSIDs()
# print(str(len(VUS_ids)) + " variants are found")


# get csv file which filters VUS_ids out with human_genes  
# df_VUS = fetchVUS_ClinVar(VUS_ids, human_genes, '../data/MM_enzyme.csv')
# print(df_VUS)

# read csv file and make dataframe
# df_VUS = readVUScsv('../data/MM_enzyme.csv')
# print(df_VUS)

# add fasta sequence to dataframe
# df_VUS = addFASTA('../data/MM_enzyme_filter_1.csv')
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

readVUS_ClinVar('../data/clinvarVariation_4.xml', '../data/MM_enzyme_short.csv', human_genes)

print('GBA' in human_genes)