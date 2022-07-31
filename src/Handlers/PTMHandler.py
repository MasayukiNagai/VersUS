from bioservices import UniProt
from collections import defaultdict
from datetime import datetime
from logging import getLogger
import xml.etree.ElementTree as ET


class PTMHandler:

    def __init__(self):
        self.uniprot = UniProt()
        self.logger = getLogger('versus_logger').getChild(__name__)

    def get_xml(self, uniprot_id):
        xml = self.uniprot.search(uniprot_id, frmt='xml')
        return xml

    def parse_xml(self, xml):
        myroot = ET.fromstring(xml)
        features = myroot[0].findall('{http://uniprot.org/uniprot}feature')
        ptm_dict = defaultdict(int)
        for feature in features:
            if feature.attrib['type'] == 'modified residue':
                pos = feature[0][0].attrib['position']
                ptm_dict[pos] += 1
        positions = ptm_dict.keys()
        return positions

    def run(self, uniprot_ids):
        self.logger.info('Start retriving PTMs')
        start = datetime.now()
        uid2ptm = {}
        for i, u_id in enumerate(uniprot_ids):
            xml = self.get_xml(u_id)
            positions = self.parse_xml(xml)
            uid2ptm[u_id] = positions
            if i % 100 == 0:
                print(f'{i}/{len(uniprot_ids)} done')
        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Retriving PTMs took {c[0]} mins {c[1]} secs')
        return uid2ptm

    def write_tsv(self, uid2ptm, outfile):
        with open(outfile, 'w') as f:
            for uid, pos_ls in uid2ptm.items():
                positions = ','.join(pos_ls)
                line = '\t'.join([uid, positions])
                f.write(line + '\n')

    def addPTM2VUSdict(self, vus_dict, gene_dict, outfile=None):
        uids = []
        for gene_id in gene_dict.keys():
            uids.append(gene_id)
        uid2ptm = self.run(uids)
        if outfile is not None:
            self.write_tsv(uid2ptm, outfile)
        for vus in vus_dict.values():
            uid = vus['uniprot_id']
            pos = vus['start']
            # end = int(vus_dict['end'])
            if pos in uid2ptm[uid]:
                vus['ptm'] = True
            else:
                vus['ptm'] = False
        return vus_dict
