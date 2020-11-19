import os
import datetime
from lxml import etree
from logging import getLogger

class BLASTHandler():

    def __init__(self, blast_path, blast_input, blast_output):
        self.blast_path = blast_path
        self.blast_input = blast_input
        self.blast_output = blast_output
        self.blast_dict = {}
        self.logger = getLogger('versus_logger').getChild(__name__)

    # make FASTA format text file from dataframe for blast search
    def make_fasta_for_blast(self, vus_dict: dict):
        # info_tup = ('NP_accession', 'gene_id', 'gene_name', 'FASTA_window')
        with open(self.blast_input, 'w') as f:
            for vus_id in range(0, len(vus_dict)):
                np_acc = vus_dict[vus_id]['NP_accession']
                gene_id = vus_dict[vus_id]['gene_id']
                gene_name = vus_dict[vus_id]['gene_name']
                seq = vus_dict[vus_id]['FASTA_window']
                seq = seq if seq != None else ''
                tup = (np_acc, gene_id, gene_name)
                header = '>' + '\t'.join(tup) + '\n'
                fasta = seq + '\n'
                f.write(header)
                f.write(fasta)


    def blast_locally(self, evalue: float=10.0):
        start = datetime.datetime.now()
        cmd1 = self.blast_path
        cmd2 = 'blastp' + ' ' \
                + '-query ' + self.blast_input + ' '\
                + '-db ' + cmd1 + 'pdbaa' + ' '\
                + '-evalue ' + str(evalue) + ' '\
                + '-outfmt ' + '5' + ' '\
                + '-out ' + self.blast_output
        cmd = os.path.join(cmd1, cmd2)
        b_cmd = os.system(cmd)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Running BLAST took {c[0]} minutes {c[1]} seconds')
        self.logger.info(cmd + ' : ran with exit code %d' %b_cmd)


    class BlastXMLHandler(object):
        def __init__(self):
            self.blast_results = {}
            self.blast_id = 0
            self.tag_stack = []
            self.info = ('pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
            self.pdb_id = ''
            self.evalue = 0
            self.hit_from = ''
            self.hit_to = ''
            self.has_got_best_pdb = False
            self.hsp_num = 0
            self.is_homo_sapiens = False
            self.has_got_species = False
            self.logger = getLogger('versus_logger').getChild(__name__)

        def start(self, tag, attrs):
            self.tag_stack.append(tag)
            if tag == 'Iteration':
                self.blast_results[self.blast_id] = {}
                
        def end(self, tag):
            self.tag_stack.pop()
            if tag == 'Hit' and not self.has_got_best_pdb:
                if (self.blast_results[self.blast_id].get('pdb_ID') == None) or \
                   ((self.evalue == self.blast_results[self.blast_id]['BLAST_evalue']) and self.is_homo_sapiens):
                    self.blast_results[self.blast_id]['pdb_ID'] = self.pdb_id
                    self.blast_results[self.blast_id]['BLAST_evalue'] = self.evalue
                    self.blast_results[self.blast_id]['hit_from'] = self.hit_from
                    self.blast_results[self.blast_id]['hit_to'] = self.hit_to
                    if self.is_homo_sapiens:
                        self.has_got_best_pdb = True
                elif self.evalue < self.blast_results[self.blast_id]['BLAST_evalue']:
                    self.has_got_best_pdb = True
                self.has_got_species = False
                
            elif tag == 'Iteration':
                for k in self.info:
                    if self.blast_results[self.blast_id].get(k) == None:
                        self.blast_results[self.blast_id][k] = None
                self.blast_id += 1
                self.has_got_best_pdb = False
                self.hsp_num = 0
                self.is_homo_sapiens = False
                if self.blast_id % 10000 == 0:
                    print(f'counter: {self.blast_id}')

        def data(self, data):
            if not self.has_got_best_pdb:
                if 'Hit_id' == self.tag_stack[-1]:
                    try:
                        pdb = data.split('pdb|')[1]
                    except: 
                        self.logger.warning(f'Cannot split {data} for pdb, id: {self.blast_id}')
                        pdb = data
                    self.pdb_id = pdb
                elif 'Hit_def' == self.tag_stack[-1] and not self.has_got_species:
                    try:
                        species = data.split('[')[1].split(']')[0]
                    except:
                        self.logger.warning(f'Cannot identify species {data}')
                        species = ''
                    if species.lower() == 'homo sapiens':
                        self.is_homo_sapiens = True
                    else:
                        self.is_homo_sapiens = False
                    self.has_got_species = True
                elif 'Hsp_num' == self.tag_stack[-1]:
                    self.hsp_num = int(data)
                elif self.hsp_num == 1:
                    if 'Hsp_evalue' == self.tag_stack[-1]:
                        self.evalue = float(data)
                    elif 'Hsp_hit-from' == self.tag_stack[-1]:
                        self.hit_from = data
                    elif 'Hsp_hit-to' == self.tag_stack[-1]:
                        self.hit_to = data
                        
        
        def close(self):
            print('Completed parsing the blast_result xml')
            print(f'Total Count: {self.blast_id}')
            print(f'Length of blast_dict: {len(self.blast_results)}')
            return self.blast_results


    def readBlastXML(self) -> dict:
        print('Start parcing')
        start = datetime.datetime.now()
        parser = etree.XMLParser(target=self.BlastXMLHandler())
        blast_dict = etree.parse(self.blast_output, parser)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        print(f'Running BlastXMLParser took {c[0]} minutes {c[1]} seconds')
        return blast_dict  


    def add_blast_results(self, vus_dict) -> dict:
        if vus_dict.keys() != self.blast_dict.keys():
            print(f'VUS_dict and Blast_dict have different keys. vus_dict: {len(vus_dict)}, blast_dict: {len(self.blast_dict)}')
            return None
        for common_id in range(0, len(vus_dict)):
            vus_dict[common_id]['pdb_ID'] = self.blast_dict[common_id]['pdb_ID']
            vus_dict[common_id]['BLAST_evalue'] = self.blast_dict[common_id]['BLAST_evalue']
            vus_dict[common_id]['hit_from'] = self.blast_dict[common_id]['hit_from']
            vus_dict[common_id]['hit_to'] = self.blast_dict[common_id]['hit_to']
        return vus_dict


    def run(self, vus_dict, evalue: float=10.0):
        self.make_fasta_for_blast(vus_dict)
        self.blast_locally(evalue)
        blast_dict = self.readBlastXML()
        vus_dict = self.add_blast_results(vus_dict)
        return vus_dict
