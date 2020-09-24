import pandas as pd
from lxml import etree
import datetime

aaMapThreeToOne = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                   'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                   'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                   'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


class VariationHandler(object):
    def __init__(self, gene_dict):
        self.vus_dict = {}
        self.gene_dict = gene_dict
        self.var_types_to_get = ('single nucleotide variant')
        self.is_var_type_to_get = False
        self.vus_id = 0
        self.clinvar_acc = ''
        self.gene_symbol = ''
        self.gene_name = ''
        self.chr = ''
        self.start_pos = ''
        self.stop_pos = ''
        self.ref = ''
        self.alt = ''
        self.np_acc = ''
        self.change = ''
        self.interpretation = ''
        self.is_first_gene_tag = True
        self.has_np_yet = False
        self.has_mut_type_yet = False
        self.is_missense = False
        self.is_uncertain = False
        self.is_conflicting = False
        self.is_not_provided = False

        # self.in_GeneList_tag = False
        # self.in_Interpretations_tag = False
        # self.in_Interpretation_tag = False
        # self.in_Description_tag = False
        # self.in_Desc_Hist_tag = False
        self.tag_stack = []
        
        self.ct_var = 0
        self.ct_missense_and_type_to_get = 0
        self.ct_uncertain_var = 0
        self.ct_conflicting_var = 0
        self.ct_not_provided_var = 0
        print('debug: start parcing')

    def start(self, tag, attrs):
        self.tag_stack.append(tag)
        if (tag == 'VariationArchive') and (attrs.get('VariationType').lower() in self.var_types_to_get):
            self.is_var_type_to_get = True  # don't forget to reset
            self.clinvar_acc = attrs.get('Accession')
        if self.is_var_type_to_get:
            # if tag == 'GeneList':
            #     self.in_GeneList_tag = True  # don't forget to reset
            if tag == 'Gene' and self.is_first_gene_tag == True:
                self.gene_symbol = attrs.get('Symbol')
                self.gene_name = attrs.get('FullName')
                self.is_first_gene_tag = False   # don't forget to reset
            # elif tag == 'SequenceLocation' and self.in_GeneList_tag == False:
            elif tag == 'SequenceLocation' and 'GeneList' not in self.tag_stack:
                if attrs.get('Assembly') == 'GRCh38':
                    self.chr = attrs.get('Chr')
                    self.start_pos = attrs.get('start')
                    self.stop_pos = attrs.get('stop')
                    self.ref = attrs.get('referenceAlleleVCF')
                    self.alt = attrs.get('alternateAlleleVCF')
            elif tag == 'ProteinExpression' and self.has_np_yet == False:
                np_acc = attrs.get('sequenceAccessionVersion')
                if (np_acc is not None) and np_acc.startswith('NP'):
                    self.np_acc = np_acc
                    self.change = attrs.get('change') 
                    self.has_np_yet = True  # don't forget to reset
            elif tag == 'MolecularConsequence' and self.has_np_yet == True and self.has_mut_type_yet == False:
                if (attrs.get('Type') is not None) and 'missense' in attrs.get('Type').lower():
                    self.is_missense = True
                    self.ct_missense_and_type_to_get += 1
                self.has_mut_type_yet = True  # don't forget to reset
        # if tag == 'Interpretations':  # if we want to skip interpretation of variants we are not interested in, we can nest the following ifs in the above if
        #     self.in_Interpretations_tag = True
        # elif tag == 'Interpretation':
        #     self.in_Interpretation_tag = True
        # elif tag == 'Description':
        #     self.in_Description_tag = True
        # elif tag == 'DescriptionHistory':
        #     self.in_Desc_Hist_tag = True

    def end(self, tag):
        self.tag_stack.pop()
        if tag == 'VariationArchive':
            clinical_significance = self.interpretation.lower()
            # print('debug: ' + clinical_significance)
            if 'uncertain' in clinical_significance:
                self.is_uncertain = True
                self.ct_uncertain_var += 1
            elif 'conflicting' in clinical_significance:
                self.is_conflicting = True
                self.ct_conflicting_var += 1
            elif 'not provided' in clinical_significance:
                self.is_not_provided = True
                self.ct_not_provided_var += 1
            if self.is_var_type_to_get \
               and (self.gene_symbol in self.gene_dict.keys())\
               and self.is_missense \
               and (self.is_uncertain or self.is_conflicting or self.is_not_provided):
                try:
                    self.change = self.change.split('p.')[1]
                    ref = aaMapThreeToOne.get(self.change[0:3])
                    alt = aaMapThreeToOne.get(self.change[len(self.change) - 3:len(self.change)])
                except:
                    ref = ''
                    alt = ''
                    # pos = int(self.change.split('_')[0][1:])
                if ref != '' and alt != '':
                    try:
                        pos = int(self.change[3:len(self.change) - 3])
                    except:
                        pos = self.change.split('_')[0][3:]
                        print(pos)
                    change_one_char = ref + str(pos) + alt
                    # register the variant
                    self.vus_dict[self.vus_id] = {}
                    self.vus_dict[self.vus_id]['ClinVar_accession'] = self.clinvar_acc
                    self.vus_dict[self.vus_id]['gene_id'] = self.gene_symbol
                    self.vus_dict[self.vus_id]['gene_name'] = self.gene_name
                    self.vus_dict[self.vus_id]['clinical_significance'] = clinical_significance
                    self.vus_dict[self.vus_id]['EC_number'] = self.gene_dict.get(self.gene_symbol)
                    self.vus_dict[self.vus_id]['missense_variation'] = change_one_char
                    self.vus_dict[self.vus_id]['NP_accession'] = self.np_acc
                    self.vus_dict[self.vus_id]['chr'] = self.chr
                    self.vus_dict[self.vus_id]['start'] = self.start_pos
                    self.vus_dict[self.vus_id]['stop'] = self.stop_pos
                    self.vus_dict[self.vus_id]['referenceAllele'] = self.ref
                    self.vus_dict[self.vus_id]['alternateAllele'] = self.alt
                    self.vus_id += 1
            # reset variables
            self.clinvar_acc = ''
            self.gene_symbol = ''
            self.gene_name = ''
            self.chr = ''
            self.start_pos = ''
            self.stop_pos = ''
            self.ref = ''
            self.alt = ''
            self.np_acc = ''
            self.change = ''
            self.interpretation = ''
            self.is_var_type_to_get = False
            self.is_first_gene_tag = True
            self.has_np_yet = False
            self.has_mut_type_yet = False
            self.is_missense = False
            self.is_uncertain = False
            self.is_conflicting = False
            self.is_not_provided = False
            self.ct_var += 1
            if self.ct_var % 10000 == 0:
                print(f'counter: {self.ct_var}')
        # elif tag == 'GeneList':
        #     self.in_GeneList_tag = False
        # elif tag == 'Interpretations':
        #     self.in_Interpretations_tag = False
        # elif tag == 'Interpretation':
        #     self.in_Interpretation_tag = False
        # elif tag == 'Description':
        #     self.in_Description_tag = False
        # elif tag == 'DescriptionHistory':
        #     self.in_Desc_Hist_tag = False
    
    def data(self, data):
        if 'Interpretations' in self.tag_stack and 'Interpretation' in self.tag_stack and 'Description' in self.tag_stack and 'DescriptionHistory' not in self.tag_stack:
            self.interpretation = data  # needs to check this part 

    def close(self):
        print(f"Total Variations: {self.ct_var}")
        print(f"Uncertain Significance: {self.ct_uncertain_var}")
        print(f"Conflicting Report: {self.ct_conflicting_var}")
        print(f"Not Provided: {self.ct_not_provided_var}")
        print(f"Missense and Type to get: {self.ct_missense_and_type_to_get}")
        print(f"VUS in the list: {len(self.vus_dict)}")
        print('debug: the file is closed')
        return self.vus_dict


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXML(xml_path, gene_dict):
    start = datetime.datetime.now()
    parser = etree.XMLParser(target=VariationHandler(gene_dict))
    vus_dict = etree.parse(xml_path, parser)
    end = datetime.datetime.now()
    time = end - start
    c = divmod(time.days * 86400 + time.seconds, 60)
    print(f'Running ClinVarXMLParser took {c[0]} minutes {c[1]} seconds')
    return vus_dict


# need to change global WFILE stuff later
class variationHandlerSpecific(object):
    def __init__(self, accession):
        self.is_accession = False
        self.accession = accession
        self.ct = 0
        print("Accession:" + self.accession)
        
    def start(self, tag, attrs):
        global WFILE
        if (tag == 'VariationArchive') and (attrs.get('Accession') == self.accession):
            self.is_accession = True
            print('The specific variation is found: ' + str(self.ct))
        if self.is_accession:
            if len(attrs.keys()) == 0:
                WFILE.write('<' + tag)
            else:   
                for i, t in enumerate(attrs.keys()):
                    if i == 0:
                        WFILE.write('<' + tag + ' ')
                    elif i != len(attrs.keys()) - 1:
                        WFILE.write(t + '="' + attrs.get(t) + '"' + ' ')
                    else:
                        WFILE.write(t + '="' + attrs.get(t) + '"')
            WFILE.write('>')
            
    def end(self, tag):
        global WFILE
        if self.is_accession:
            WFILE.write('</' + tag + '>')
        if tag == 'VariationArchive':
            if self.ct % 10000 == 0:
                print(self.ct)
            self.ct += 1
        if self.is_accession and tag == 'VariationArchive':
            self.is_accession = False
            WFILE.close()
            print('The subnode file is completed')
            
    def data(self, data):
        global WFILE
        if data is not None:
            if self.is_accession and data != "":
                WFILE.write(data)
            
    def close(self):
        print('The xml file is closed')


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXMLSpecific(input_path, accession):
    print('Start parcing')
    parser = etree.XMLParser(target=variationHandlerSpecific(accession))
    etree.parse(input_path, parser)


class BlastHandler(object):
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


def readBlastXML(input_path):
    print('Start parcing')
    start = datetime.datetime.now()
    parser = etree.XMLParser(target=BlastHandler())
    blast_dict = etree.parse(input_path, parser)
    end = datetime.datetime.now()
    time = end - start
    c = divmod(time.days * 86400 + time.seconds, 60)
    print(f'Running BlastXMLParser took {c[0]} minutes {c[1]} seconds')
    return blast_dict    
