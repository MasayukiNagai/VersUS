import os
import csv
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "mnaffairs.intl@gmail.com"

class SeqHandler:

    def __init__(self, gene_file, seq_dir):
        self.gene_file = gene_file
        self.genes_dict = {}
        self.seq_dir = seq_dir
        self.seq_dict = {}

    
    # read csv file of gene names of human enzymes
    # return dict{gene_names: EC number}
    def readHumanGenesEC(self):
        dup = set()
        with open(self.gene_file, 'r') as f:
            f.readline()
            reader = csv.reader(f)
            for row in reader:
                key = row[0]
                if key in self.genes_dict.keys():
                    dup.add(key)
                else:
                    self.genes_dict[key] = row[1]
        print(f'{len(dup)} genes are duplicated: ', dup)  # prints, if any, genes duplicated in the file
        return self.genes_dict


    # makes dictionary of fasta sequences and np number 
    # returns dict{np_num: fasta sequence}
    def make_seq_dict(self):
        ct_np = 0
        for root, d_names, file_names in os.walk(self.seq_dir):
            for filename in file_names:
                fname = os.path.join(root, filename)
                with open(fname, 'r') as f:
                    print(f'opened a fasta file: {fname}')
                    np_num = ''
                    for line in f:
                        if '>' in line:
                            ct_np += 1
                            try:
                                np_num = line.split()[0].split('>')[1]
                                self.seq_dict[np_num] = ''
                            except:
                                print(f'Error: {line}')
                                return
                            if line[0] != '>':
                                print(f'Not start with >: {line} in {fname}')
                        else:
                            seq_dict[np_num] += line.strip()
        print(f'The number of items in the seq dict: {len(self.seq_dict)}')
        print(f'The number of np found in the files: {ct_np}')
        return self.seq_dict


    # fetches fasta sequences for varinants whose sequences aren't in the imported files
    # returns dict{np_num: fasta sequence}
    def fetch_seq(self, ls_np):
        for np_num in ls_np:
            handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
            seq_record = SeqIO.read(handle, 'fasta')
            sequence = seq_record.seq
            self.seq_dict[np_num] = sequence
        return self.seq_dict


    # crops fasta sequence
    # returns the sequnece cropped in a specified range
    def crop_seq(self, seq: str, pos: int, ref: str, seq_range: int) -> str:
        if pos-1 < len(seq) and seq[pos-1].upper() == ref.upper():
            if pos-1 - seq_range <= 0:
                proteinSeq = seq[0:2*seq_range+1]
            elif pos-1 + seq_range > len(seq) - 1:
                proteinSeq = seq[len(seq)-1-2*seq_range if len(seq)-1-2*seq_range >= 0 else 0:len(seq)]
            else:
                proteinSeq = seq[pos-1-seq_range:pos+seq_range]
            return proteinSeq
        else:
            print('Error: crop_seq')
            return None
    
    
    # add fasta sequence to vus dict
    # return vus dict and unfound_seq set
    def add_seq_to_dict(self, vus_dict: dict, seq_range: int=12):
        unfound_seq = set()
        seq_list = []
        for vus_id in vus_dict.keys():
            if vus_dict[vus_id].get('FASTA_window') != None:
                continue
            mutation = vus_dict[vus_id]['missense_variation']
            ref = mutation[0]
            pos = int(mutation[1:-1])
            np_num = vus_dict[vus_id]['NP_accession']
            try:
                seq = self.seq_dict[np_num]
                # seq_cropped = ''
                seq_cropped = self.crop_seq(seq, pos, ref, seq_range)
            except:
                seq_cropped = None
                unfound_seq.add(np_num)
            # if(vus_id==1):
            #     print(mutation, ref, pos, np_num, seq)
            vus_dict[vus_id]['FASTA_window'] = seq_cropped
        return vus_dict, unfound_seq