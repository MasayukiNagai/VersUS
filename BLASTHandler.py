import os
import datetime
from lxml import etree

class BLASTHandler():

    def __init__(self, blast_path, blast_input, blast_output):
        self.blast_path = blast_path
        self.blast_input = blast_input
        self.blast_output = blast_output
        self.blast_dict = {}


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
        print(f'Running BLAST took {c[0]} minutes {c[1]} seconds')
        print(cmd + ' : ran with exit code %d' %b_cmd)


    def add_blast_results(self, vus_dict):
        if vus_dict.keys() != self.blast_dict.keys():
            print(f'VUS_dict and Blast_dict have different keys. vus_dict: {len(vus_dict)}, blast_dict: {len(self.blast_dict)}')
            return None
        for common_id in range(0, len(vus_dict)):
            vus_dict[common_id]['pdb_ID'] = self.blast_dict[common_id]['pdb_ID']
            vus_dict[common_id]['BLAST_evalue'] = self.blast_dict[common_id]['BLAST_evalue']
            vus_dict[common_id]['hit_from'] = self.blast_dict[common_id]['hit_from']
            vus_dict[common_id]['hit_to'] = self.blast_dict[common_id]['hit_to']
        return vus_dict


    class BlastXMLHandler(object):
        def __init__(self):
            self.blast_results = {}
            self.blast_id = 0
            self.tag_stack = []
            self.info = ('pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
            # self.dictlist = {'pdb_ID': [], 'BLAST_evalue': [], 'hit_from': [], 'hit_to': []}
            # self.is_hit_id = False
            # self.is_evalue = False
            # self.is_hit_from = False
            # self.is_hit_to = False
            self.ct_hit = 0

        def start(self, tag, attrs):
            self.tag_stack.append(tag)
            if tag == 'Iteration':
                self.blast_results[self.blast_id] = {}
            elif tag == 'Hit':  # there are a few <Hsp> tags within a <Hit> tag
                self.ct_hit += 1
            # elif tag == 'Hit_id':
            #     self.is_hit_id = True
            # elif tag == 'Hsp_evalue':
            #     self.is_evalue = True
            # elif tag == 'Hsp_hit-from':
            #     self.is_hit_from = True
            # elif tag == 'Hsp_hit-to':
            #     self.is_hit_to = True
                
        def end(self, tag):
            self.tag_stack.pop()
            if tag == 'Iteration':
                # if self.ct_hit == 0:  # when there is no hit, add None to each list
                #     self.blast_results[self.blast_id] = {}
                #     self.blast_results[self.blast_id]['pdb_ID'] = None
                #     self.blast_results[self.blast_id]['BLAST_evalue'] = None
                #     self.blast_results[self.blast_id]['hit_from'] = None
                #     self.blast_results[self.blast_id]['hit_to'] = None
                for k in self.info:
                    if self.blast_results[self.blast_id].get(k) == None:
                        self.blast_results[self.blast_id][k] = None
                self.blast_id += 1
                self.ct_hit = 0
                if self.blast_id % 10000 == 0:
                    print(f'counter: {self.blast_id}')
            elif tag == 'Hsp':
                self.ct_hit += 1
            # elif tag == 'Hit_id':
            #     self.is_hit_id = False
            # elif tag == 'Hsp_evalue':
            #     self.is_evalue = False
            # elif tag == 'Hsp_hit-from':
            #     self.is_hit_from = False
            # elif tag == 'Hsp_hit-to':
            #     self.is_hit_to = False

        def data(self, data):
            if self.ct_hit == 1:
                if 'Hit_id' in self.tag_stack:
                    try:
                        pdb = data.split('pdb|')[1]
                    except: 
                        print(f'Cannot split {data} for pdb, id: {self.blast_id}')
                        pdb = data
                    self.blast_results[self.blast_id]['pdb_ID'] = pdb
                elif 'Hsp_evalue' in self.tag_stack:
                    self.blast_results[self.blast_id]['BLAST_evalue'] = data
                elif 'Hsp_hit-from' in self.tag_stack:
                    self.blast_results[self.blast_id]['hit_from'] = data
                elif 'Hsp_hit-to' in self.tag_stack:
                    self.blast_results[self.blast_id]['hit_to'] = data
        
        def close(self):
            print('Completed parsing the blast_result xml')
            print(f'Total Count: {self.blast_id}')
            print(f'Length of blast_dict: {len(self.blast_results)}')
            return self.blast_results


    def readBlastXML(self):
        print('Start parcing')
        start = datetime.datetime.now()
        parser = etree.XMLParser(target=BlastXMLHandler())
        blast_dict = etree.parse(self.blast_output, parser)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        print(f'Running BlastXMLParser took {c[0]} minutes {c[1]} seconds')
        return blast_dict  

