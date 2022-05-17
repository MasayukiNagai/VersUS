from this import d
from lxml import etree
import gzip
from collections import defaultdict
from datetime import datetime
from logging import getLogger

aaMapThreeToOne = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                   'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                   'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                   'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

class ClinVarHandler:

    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml
        self.logger = getLogger('versus_logger').getChild(__name__)


    class VariationHandler(object):
        def __init__(self, gene_set):
            self.logger = getLogger('versus_logger').getChild(__name__)
            self.vus_dict = {}
            self.gene_set = gene_set
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

            self.tag_stack = []

            self.ct_var = 0
            self.ct_missense_and_type_to_get = 0
            self.ct_uncertain_var = 0
            self.ct_conflicting_var = 0
            self.ct_not_provided_var = 0

            self.var_types = defaultdict(int)
            self.clinical_significances = defaultdict(int)
            self.change_types = defaultdict(int)
            self.is_enzyme = defaultdict(int)

        def start(self, tag, attrs):
            self.tag_stack.append(tag)
            if (tag == 'VariationArchive'):
                self.var_types[attrs.get('VariationType')] += 1
                if (attrs.get('VariationType').lower() in self.var_types_to_get):
                    self.is_var_type_to_get = True
                    # self.clinvar_acc = attrs.get('Accession')
                    self.clinvar_acc = attrs.get('VariationID')
            if self.is_var_type_to_get:
                if tag == 'Gene' and self.is_first_gene_tag == True:
                    self.gene_symbol = attrs.get('Symbol')
                    self.gene_name = attrs.get('FullName')
                    self.is_first_gene_tag = False
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
                        self.has_np_yet = True
                elif tag == 'MolecularConsequence' and self.has_np_yet == True and self.has_mut_type_yet == False:
                    if (attrs.get('Type') is not None) and 'missense' in attrs.get('Type').lower():
                        self.is_missense = True
                        self.ct_missense_and_type_to_get += 1
                    self.has_mut_type_yet = True

        def end(self, tag):
            self.tag_stack.pop()
            if tag == 'VariationArchive':
                self.clinical_significances[self.interpretation] += 1
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
                if self.gene_symbol in self.gene_set:
                    self.is_enzyme['yes'] += 1
                else:
                    self.is_enzyme['no'] += 1
                if self.is_var_type_to_get \
                   and (self.gene_symbol in self.gene_set)\
                   and self.is_missense\
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
                if self.ct_var % 100000 == 0:
                    print(f'{self.ct_var} variations have been processed')

        def data(self, data):
            if 'Interpretations' in self.tag_stack and 'Interpretation' in self.tag_stack and 'Description' in self.tag_stack and 'DescriptionHistory' not in self.tag_stack:
                self.interpretation = data  # needs to check this part

        def close(self):
            self.logger.info(f'Finish parsing ClinvarVariaiton\n\
                               Total Variations: {self.ct_var}\n\
                               Uncertain Significance: {self.ct_uncertain_var}\n\
                               Conflicting Report: {self.ct_conflicting_var}\n\
                               Not Provided: {self.ct_not_provided_var}\n\
                               Missense with specified type(s): {self.ct_missense_and_type_to_get}\n\
                               VUS in the list: {len(self.vus_dict)}')
            self.logger.info(f'Clinical Significances: {self.clinical_significances}\n\
                               Variable Types: {self.var_types}\n\
                               Enzyme Ratio: {self.is_enzyme}')
            return self.vus_dict


    # read xml file of variations from ClinVar
    # return dataframe and write to a csv file
    def readClinVarVariationsXML(self, gene_set):
        self.logger.info('Start parcing ClinvarVariationsRelease')
        start = datetime.now()
        parser = etree.XMLParser(target=self.VariationHandler(gene_set))
        with gzip.open(self.clinvar_xml, 'rt') as handle:
            vus_dict = etree.parse(handle, parser)
        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Running ClinVarXMLParser took {c[0]} minutes {c[1]} seconds')
        return vus_dict


    # # need to change global WFILE stuff later
    # class VariationHandlerSpecific(object):
    #     def __init__(self, accession):
    #         self.is_accession = False
    #         self.accession = accession
    #         self.ct = 0
    #         print("Accession:" + self.accession)

    #     def start(self, tag, attrs):
    #         global WFILE
    #         if (tag == 'VariationArchive') and (attrs.get('Accession') == self.accession):
    #             self.is_accession = True
    #             print('The specific variation is found: ' + str(self.ct))
    #         if self.is_accession:
    #             if len(attrs.keys()) == 0:
    #                 WFILE.write('<' + tag)
    #             else:
    #                 for i, t in enumerate(attrs.keys()):
    #                     if i == 0:
    #                         WFILE.write('<' + tag + ' ')
    #                     elif i != len(attrs.keys()) - 1:
    #                         WFILE.write(t + '="' + attrs.get(t) + '"' + ' ')
    #                     else:
    #                         WFILE.write(t + '="' + attrs.get(t) + '"')
    #             WFILE.write('>')

    #     def end(self, tag):
    #         global WFILE
    #         if self.is_accession:
    #             WFILE.write('</' + tag + '>')
    #         if tag == 'VariationArchive':
    #             if self.ct % 10000 == 0:
    #                 print(self.ct)
    #             self.ct += 1
    #         if self.is_accession and tag == 'VariationArchive':
    #             self.is_accession = False
    #             WFILE.close()
    #             print('The subnode file is completed')

    #     def data(self, data):
    #         global WFILE
    #         if data is not None:
    #             if self.is_accession and data != "":
    #                 WFILE.write(data)

    #     def close(self):
    #         print('The xml file is closed')


    # # read xml file of variations from ClinVar
    # # return dataframe and write to a csv file
    # def readClinVarVariationsXMLSpecific(self, accession):
    #     self.logger.info('Start parcing ClinVarVariations Relase (Specific)')
    #     parser = etree.XMLParser(target=self.VariationHandlerSpecific(accession))
    #     etree.parse(self.clinvar_xml, parser)


    def run(self, gene_set):
        vus_dict = self.readClinVarVariationsXML(gene_set)
        self.logger.info('Finish processing the ClinVarVariationRelease')
        return vus_dict
