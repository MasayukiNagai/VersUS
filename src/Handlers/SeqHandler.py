import os
import gzip
from Bio import Entrez
from Bio import SeqIO
from logging import getLogger



class SeqHandler:

    def __init__(self):
        self.logger = getLogger('versus_logger').getChild(__name__)

    def setup(self, email, api_key):
        Entrez.email = email
        self.api_key = api_key

    # read tsv file of entry(uniprot_id) - gene id - ec number
    # return dict {gene_id: 'ec': ec number, 'uniprot_id': uniprot_id}
    def readUniprot_GeneId_EC(self, genefile):
        dup = set()
        genes_dict = {}
        with open(genefile, 'r') as f:
            f.readline()  # skips the header
            for line in f:
                data = line.rstrip().split('\t')
                uniprot_id = data[0]
                gene_id = data[1]
                ec = data[2]
                if gene_id in genes_dict.keys():
                    dup.add(gene_id)
                else:
                    genes_dict[gene_id] = {'ec': ec, 'uniprot_id': uniprot_id}
        self.logger.info(f'{len(dup)} genes are duplicated: {dup}')
        return genes_dict

    def add_uniprotId_EC(self, vus_dict, genes_dict):
        uids = set()
        for vus in vus_dict.values():
            gene_id = vus['gene_id']
            uniprot_id = genes_dict[gene_id]['uniprot_id']
            ec = genes_dict[gene_id]['ec']
            vus['uniprot_id'] = uniprot_id
            vus['EC_number'] = ec
            uids.add(uniprot_id)
        return vus_dict, uids

    # makes dictionary of fasta sequences and np number
    # returns dict{np_num: fasta sequence}
    def make_seq_dict(self, proteomesdir, seq_dict={}):
        self.logger.debug(f'Open sequence files')
        ct_np = 0
        for root, _, file_names in os.walk(proteomesdir):
            for filename in file_names:
                fname = os.path.join(root, filename)
                with gzip.open(fname, 'rt') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        seq_dict[record.id] = str(record.seq)
        self.logger.debug(f'Finish storing {len(seq_dict)} sequences')
        # print(f'The number of np found in the files: {ct_np}')
        return seq_dict

    # fetches fasta sequences for varinants whose sequences aren't in the imported files
    # returns dict{np_num: fasta sequence}
    def fetch_seq(self, nps: list):
        seq_dict = {}
        for np_num in nps:
            handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta',
                                   retmode='text', api_key=self.api_key)
            seq_record = SeqIO.read(handle, 'fasta')
            seq_dict[np_num] = str(seq_record.seq)
        return seq_dict

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
    def add_seq_to_dict(self, vus_dict: dict, seq_dict: dict,  seq_range: int=12):
        unfound_seq = set()
        seq_list = []
        for vus in vus_dict.values():
            if vus.get('FASTA_window') is not None:
                continue
            mutation = vus['missense_variation']
            ref = mutation[0]
            pos = int(mutation[1:-1])
            np_num = vus['NP_accession']
            try:
                seq = seq_dict[np_num]
                seq_cropped = self.crop_seq(seq, pos, ref, seq_range)
            except:
                seq_cropped = None
                unfound_seq.add(np_num)
            vus['FASTA_window'] = seq_cropped
        return vus_dict, unfound_seq

    def get_seq(self, vus_dict, proteomesdir, window_size: int=12):
        self.logger.info('Add sequences to vus_dict')
        seq_dict = self.make_seq_dict(proteomesdir)
        vus_dict, unfonund_np_ls = self.add_seq_to_dict(
            vus_dict, seq_dict, window_size)
        seq_fetched_dict = self.fetch_seq(unfonund_np_ls)
        vus_dict, _ = self.add_seq_to_dict(
            vus_dict, seq_fetched_dict, window_size)
        self.logger.info('Finish adding sequences to vus_dict')
        return vus_dict
