import pandas as pd
import numpy as np
import csv
import os
from Bio import Entrez
from Bio import SeqIO
# from Bio.Blast import NCBIWWW
import time
import datetime

Entrez.email = "mnaffairs.intl@gmail.com"

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