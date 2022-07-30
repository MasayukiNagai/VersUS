from bioservices import UniProt
from collections import defaultdict
import xml.etree.ElementTree as ET


class PTMHandler:

    def __init__(self):
        self.uniprot = UniProt()

    def get_xml(self, uniprot_id):
        xml = self.uniprot.serach(uniprot_id, frmt='xml')
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
        uid2ptm = {}
        for u_id in uniprot_ids:
            xml = self.get_xml(u_id)
            positions = self.parse_xml(xml)
            uid2ptm[u_id] = positions
        return uid2ptm

    def addPTM2VUSdict(self, vus_dict, gene_dict):
        uids = []
        for gene_id in gene_dict.keys():
            uids.append(gene_id)
        uid2ptm = self.run(uids)
        for vus in vus_dict.values():
            uid = vus['uniprot_id']
            pos = vus['start']
            # end = int(vus_dict['end'])
            if pos in uid2ptm[uid]:
                vus['ptm'] = True
            else:
                vus['ptm'] = False
        return vus_dict