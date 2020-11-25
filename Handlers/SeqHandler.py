import os
import csv
from Bio import Entrez
from Bio import SeqIO
from logging import getLogger

Entrez.email = "mnaffairs.intl@gmail.com"

class SeqHandler:

    def __init__(self, gene_file, proteomes_dir):
        self.gene_file = gene_file
        self.genes_dict = {}
        self.proteomes_dir = proteomes_dir
        self.seq_dict = {}
        self.logger = getLogger('versus_logger').getChild(__name__)

    
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
        self.logger.info(f'{len(dup)} genes are duplicated: {dup}')
        return self.genes_dict


    # makes dictionary of fasta sequences and np number 
    # returns dict{np_num: fasta sequence}
    def make_seq_dict(self):
        self.logger.debug(f'Open sequence files')
        ct_np = 0
        for root, d_names, file_names in os.walk(self.proteomes_dir):
            for filename in file_names:
                fname = os.path.join(root, filename)
                for record in SeqIO.parse(fname, 'fasta'):
                    self.seq_dict[record.id] = str(record.seq)
        self.logger.debug(f'Finish storing {len(self.seq_dict)} sequences')
        # print(f'The number of np found in the files: {ct_np}')
        return self.seq_dict


    # fetches fasta sequences for varinants whose sequences aren't in the imported files
    # returns dict{np_num: fasta sequence}
    def fetch_seq(self, ls_np):
        for np_num in ls_np:
            handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
            seq_record = SeqIO.read(handle, 'fasta')
            self.seq_dict[np_num] = str(seq_record.seq)
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
            self.logger.warning('Failed to crop a sequence')
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
    
    
    def get_seq(self, vus_dict, window_size: int=12):
        self.logger.info('Add sequences to vus_dict')
        self.make_seq_dict()
        vus_dict, unfonund_np_ls = self.add_seq_to_dict(vus_dict, window_size)
        self.fetch_seq(unfonund_np_ls)
        vus_dict, unfonund_np_ls = self.add_seq_to_dict(vus_dict, window_size)
        self.logger.info('Finish adding sequences to vus_dict')
        return vus_dict
