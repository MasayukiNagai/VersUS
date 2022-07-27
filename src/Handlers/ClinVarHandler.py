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
            self.refAllele = ''
            self.altAllele = ''
            self.hgvs_ls = []
            self.is_MANESelect = False
            self.has_ProteinExpression = False
            self.has_MANEandMissense = False
            # self.np_acc = ''
            # self.change = ''
            self.interpretation = ''
            self.is_first_gene_tag = True
            # self.has_np_yet = False
            # self.has_mut_type_yet = False
            # self.is_missense = False
            # self.is_uncertain = False
            # self.is_conflicting = False
            # self.is_not_provided = False

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
                        self.refAllele = attrs.get('referenceAlleleVCF')
                        self.altAllele = attrs.get('alternateAlleleVCF')
                elif not self.has_MANEandMissense:
                    if tag == 'NucleotideExpression':
                        if attrs.get('MANESelect'):
                            self.is_MANESelect = True
                    elif tag == 'ProteinExpression':
                        np_acc = attrs.get('sequenceAccessionVersion')
                        if np_acc and np_acc.startswith('NP'):
                            hgvs = {'NP': np_acc, 'change': attrs.get('change')}
                            self.hgvs_ls.append(hgvs)
                            self.has_ProteinExpression = True
                    elif tag == 'MolecularConsequence':
                        if self.has_ProteinExpression:
                            is_missense = 'missense' in attrs.get('Type').lower()
                            hgvs = self.hgvs_ls[-1]
                            hgvs['missense'] = is_missense or hgvs.get('missense')
                            if self.is_MANESelect and is_missense:
                                self.has_MANEandMissense = True
                # elif tag == 'MolecularConsequence' and self.has_np_yet == True and self.has_mut_type_yet == False:
                #     if attrs.get('Type') is not None:
                #         self.change_types[attrs.get('Type')] += 1
                #     else:
                #         self.change_types['None'] += 1
                #     if (attrs.get('Type') is not None) and 'missense' in attrs.get('Type').lower():
                #         self.is_missense = True
                #         self.ct_missense_and_type_to_get += 1
                #     self.has_mut_type_yet = True

        def end(self, tag):
            self.tag_stack.pop()
            if tag == 'HGVS':
                self.is_MANESelect = False
                self.has_ProteinExpression = False
            elif tag == 'VariationArchive':
                self.clinical_significances[self.interpretation] += 1
                # Check clinical significance
                clinical_significance = self.interpretation.lower()
                is_uncertain, is_conflicting, is_not_provided = [False] * 3
                if 'uncertain' in clinical_significance:
                    is_uncertain = True
                    self.ct_uncertain_var += 1
                elif 'conflicting' in clinical_significance:
                    is_conflicting = True
                    self.ct_conflicting_var += 1
                elif 'not provided' in clinical_significance:
                    is_not_provided = True
                    self.ct_not_provided_var += 1
                if self.gene_symbol in self.gene_set:
                    self.is_enzyme['yes'] += 1
                else:
                    self.is_enzyme['no'] += 1
                # Check missense
                is_missense = False
                if self.has_MANEandMissense:
                    hgvs = self.hgvs_ls[-1]
                    is_missense = True
                else:
                    for v in self.hgvs_ls:
                        if v.get('missense'):
                            hgvs = v
                            is_missense = True
                if is_missense:
                    self.ct_missense_and_type_to_get += 1
                # Check the whole criteria
                if self.is_var_type_to_get \
                   and (self.gene_symbol in self.gene_set)\
                   and is_missense\
                   and (is_uncertain or is_conflicting or is_not_provided):
                    try:
                        change = hgvs['change'].split('p.')[1]
                        ref = aaMapThreeToOne.get(change[0:3])
                        alt = aaMapThreeToOne.get(change[len(change) - 3:len(change)])
                    except:
                        ref = None
                        alt = None
                        self.logger.debug(f'Failed to parse "{change}"')
                    if ref and alt:
                        try:
                            pos = int(change[3:len(change) - 3])
                        except:
                            pos = change.split('_')[0][3:]
                            print(pos)
                        change_one_char = ref + str(pos) + alt
                        # Register the variant
                        self.vus_dict[self.vus_id] = {}
                        self.vus_dict[self.vus_id]['ClinVar_accession'] = self.clinvar_acc
                        self.vus_dict[self.vus_id]['gene_id'] = self.gene_symbol
                        self.vus_dict[self.vus_id]['gene_name'] = self.gene_name
                        self.vus_dict[self.vus_id]['clinical_significance'] = clinical_significance
                        self.vus_dict[self.vus_id]['missense_variation'] = change_one_char
                        self.vus_dict[self.vus_id]['NP_accession'] = hgvs['NP']
                        self.vus_dict[self.vus_id]['chr'] = self.chr
                        self.vus_dict[self.vus_id]['start'] = self.start_pos
                        self.vus_dict[self.vus_id]['stop'] = self.stop_pos
                        self.vus_dict[self.vus_id]['referenceAllele'] = self.refAllele
                        self.vus_dict[self.vus_id]['alternateAllele'] = self.altAllele
                        self.vus_id += 1
                # Reset variables
                self.clinvar_acc = ''
                self.gene_symbol = ''
                self.gene_name = ''
                self.chr = ''
                self.start_pos = ''
                self.stop_pos = ''
                self.refAllele = ''
                self.altAllele = ''
                # self.np_acc = ''
                # self.change = ''
                self.hgvs_ls = []
                self.has_MANEandMissense = False
                self.interpretation = ''
                self.is_var_type_to_get = False
                self.is_first_gene_tag = True
                # self.has_np_yet = False
                # self.has_mut_type_yet = False
                # self.is_missense = False
                # self.is_uncertain = False
                # self.is_conflicting = False
                # self.is_not_provided = False
                self.ct_var += 1
                if self.ct_var % 100000 == 0:
                    print(f'{self.ct_var} variations have been processed')

        def data(self, data):
            if 'Interpretations' in self.tag_stack and 'Interpretation' in self.tag_stack and 'Description' in self.tag_stack and 'DescriptionHistory' not in self.tag_stack:
                self.interpretation = data  # needs to check this part

        def close(self):
            self.logger.info(
                (f'Finish parsing ClinvarVariaiton\n'
                 f'Total Variations: {self.ct_var}\n'
                 f'Uncertain Significance: {self.ct_uncertain_var}\n'
                 f'Conflicting Report: {self.ct_conflicting_var}\n'
                 f'Not Provided: {self.ct_not_provided_var}\n'
                 f'Missense with specified type(s): {self.ct_missense_and_type_to_get}\n'
                 f'VUS in the list: {len(self.vus_dict)}'))
            self.logger.info(
                (f'Clinical Significances: {self.clinical_significances}\n'
                 f'Variable Types: {self.var_types}\n'
                 f'Enzyme Ratio: {self.is_enzyme}\n'
                 f'Change Types: {self.change_types}'))
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


    def run(self, gene_set):
        vus_dict = self.readClinVarVariationsXML(gene_set)
        self.logger.info('Finish processing the ClinVarVariationRelease')
        return vus_dict
