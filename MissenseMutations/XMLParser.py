import pandas as pd
from lxml import etree

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
        self.clinvar_acc = ''
        self.gene_symbol = ''
        self.gene_name = ''
        self.chr = ''
        self.start_pos = ''
        self.stop_pos = ''
        self.ref = ''
        self.alt = ''
        self.is_first_gene_tag = True
        self.has_np_yet = False
        self.has_mut_type_yet = False
        self.is_missense = False
        self.is_uncertain = False
        self.is_conflicting = False
        self.is_not_provided = False

        self.in_GeneList_tag = False
        self.in_Interpretations_tag = False
        self.in_Interpretation_tag = False
        self.in_Description_tag = False
        self.in_Desc_Hist_tag = False
        
        self.ct_var = 0
        self.ct_missense_var = 0
        self.ct_uncertain_var = 0
        self.ct_conflicting_var = 0
        self.ct_not_provided_var = 0
        print('debug: start parcing')

    def start(self, tag, attrs):
        if (tag == 'VariationArchive') and (attrs.get('VariationType').lower() in self.var_types_to_pick):
            self.is_var_type_to_get = True  # don't forget to reset
            self.clinvar_acc = attrs.get('Accession')
        if self.is_var_type_to_get:
            if tag == 'GeneList':
                self.in_GeneList_tag = True  # don't forget to reset
            elif tag == 'Gene' and self.is_first_gene_tag == True:
                self.gene_symbol = attrs.get('Symbol')
                self.gene_name = attrs.get('FullName')
                self.is_first_gene_tag == False   # don't forget to reset
            elif tag == 'SequenceLocation' and self.in_GeneList_tag == False:
                if attrs.get('Assembly') == 'GRCh38':
                    self.chr = attrs.get('Chr')
                    self.start_pos = attrs.get('start')
                    self.stop_pos = attrs.get('stop')
                    self.ref = attrs.get('referenceAlleleVCF')
                    self.alt = attrs.get('alternateAlleleVCF')
            elif tag == 'ProteinExpression' and self.has_np_yet == False:
                np_acc = attrs.get('sequenceAccessionVersion')
                if (np_acc is not None) and self.np_acc.startswith('NP'):
                    self.np_acc = np_acc
                    self.change = attrs.get('change') 
                    self.has_np_yet = True  # don't forget to reset
            elif tag == 'MolecularConsequence' and self.has_np_yet == True and self.has_mut_type_yet == False:
                if (attrs.get('Type') is not None) and 'missense' in attrs.get('Type').lower():
                    self.is_missense = True
                    self.ct_missense_var += 1
                self.has_mut_type_yet = True  # don't forget to reset
            elif tag == 'Interpretations':
                self.in_Interpretations_tag = True
            elif tag == 'Interpretation':
                self.in_Interpretation_tag = True
            elif tag == 'Description':
                self.in_Description_tag = True
            elif tag == 'DescriptionHistory':
                self.in_Desc_Hist_tag = True

    def end(self, tag):
        if tag == 'VariationArchive' and self.is_var_type_to_get:
            clinical_significance = self.interpretation.lower()
            if 'uncertain' in clinical_significance:
                self.is_uncertain = True
                self.ct_uncertain_var += 1
            elif 'conflicting' in clinical_significance:
                self.is_conflicting = True
                self.ct_conflicting_var += 1
            elif 'not provided' in clinical_significance:
                self.is_not_provided = True
                self.ct_not_provided += 1
            if (self.gene_symbol in self.gene_dict.keys()) and self.is_missense and (self.is_uncertain or self.is_conflicting or self.is_not_provided):
                try:
                    self.change = self.change.split('p.')[1]
                    ref = aaMapThreeToOne.get(self.change[0:3])
                    alt = aaMapThreeToOne.get(self.change[len(self.change) - 3:len(self.change)])
                except:
                    ref = None
                    alt = None
                if ref and alt:
                    pos = self.change[3:len(self.change) - 3]
                    change_one_char = ref + pos + alt
                    self.vus_dict[self.clinvar_acc] = {}
                    self.vus_dict[self.clinvar_acc]['gene_symbol'] = self.gene_symbol
                    self.vus_dict[self.clinvar_acc]['gene_name'] = self.gene_name
                    self.vus_dict[self.clinvar_acc]['clinical_significance'] = clinical_significance
                    self.vus_dict[self.clinvar_acc]['EC_number'] = self.gene_dict.get(self.gene_symbol)
                    self.vus_dict[self.clinvar_acc]['protein_change'] = change_one_char
                    self.vus_dict[self.clinvar_acc]['chr'] = self.chr
                    self.vus_dict[self.clinvar_acc]['start'] = self.start_pos
                    self.vus_dict[self.clinvar_acc]['stop'] = self.stop_pos
                    self.vus_dict[self.clinvar_acc]['referenceAllele'] = self.ref
                    self.vus_dict[self.clinvar_acc]['alternateAllele'] = self.alt

        if tag == 'VariationArchive':
            self.is_var_type_to_get = False
            self.is_first_gene_tag = True
            self.has_np_yet = False
            self.has_mut_type_yet = False
            self.is_missense = False
            self.is_uncertain = False
            self.is_conflicting = False
            self.is_not_provided = False
            self.ct_var += 1
            if self.ct % 10000 == 0:
                print(f'counter: {self.ct_var}')
        elif tag == 'GeneList':
            self.in_GeneList_tag = False
        elif tag == 'Interpretations':
            self.in_Interpretations_tag = False
        elif tag == 'Interpretation':
            self.in_Interpretation_tag = False
        elif tag == 'Description':
            self.in_Description_tag = False
        elif tag == 'DescriptionHistory':
            self.in_Desc_Hist_tag = False
    
    def data(self, data):
        if self.in_Interpretations_tag and self.in_Interpretation_tag and self.in_Description_tag and (not self.in_Desc_Hist_tag):
            self.interpretation = data  # needs to check this part 

    def close(self):
        print(f"Variations: {self.ct_var}")
        print(f"Uncertain Significance: {self.ct_uncertain_var}")
        print(f"Conflicting Report: {self.ct_conflicting_var}")
        print(f"Not Provided: {self.ct_not_provided_var}")
        print(f"Missense: {self.ct_missense_var}")
        print(f"Mutations in the list: {len(self.vus_dict)}")
        print('debug: the file is closed')
        return self.vus_dict


class variationHandler(object):
    def __init__(self, genes_dict):
        self.dictlist = {'gene_ID':[], 'gene_name':[], 'clinical_significance': [], 'EC_number': [], 'missense_variation': [], 'NP_accession': [], 'ClinVar_accession':[], 'Chr': [], 'start':[], 'stop':[], 'referenceAllele':[], 'alternateAllele':[]}
        self.genes_dict = genes_dict
        self.unnecessary_types = ('inversion', 'copy number gain', 'tandem duplication', 'microsatellite', 'copy number loss', 'distinct chromosomes', 'fusion', 'complex', 'duplication', 'translocation')
        self.is_type = False
        self.gene_ID = ""
        self.gene_name = ""
        self.clinvar_acc = ""
        self.np_acc = ""
        self.change = ""
        self.chr = ""
        self.start_num = ""
        self.stop_num = ""
        self.referenceAllele = ""
        self.alternateAllele = ""
        self.is_GeneList = False
        self.ct_gene = 0
        self.check_grch = False
        self.ct_np = 0
        self.ct_mc = 0  # counter for Molecular Consequence tag
        self.is_haplotype = False  # check if a variation is haplotype or not
        self.is_genotype = False  # check if a variation is genotype or not
        self.is_missense = False
        self.is_conflicting = False
        self.is_not_provided = False
        self.is_interpretations = False
        self.is_interpretation = False
        self.is_description = False
        self.is_desc_hist = False
        self.intpn = []
        self.ct = 0
        self.ct_missense = 0
        self.ct_uncertain = 0
        self.ct_conflicting = 0
        self.ct_not_provided = 0
        
    def start(self, tag, attrs):
        if (tag == 'VariationArchive') and (attrs.get('VariationType').lower() not in self.unnecessary_types):
            self.is_type = True
            if attrs.get('VariationType').lower() == 'haplotype':
                self.is_haplotype = True
            if attrs.get('VariationType').lower() == 'compoundheterozygote':
                self.is_genotype = True
            self.clinvar_acc = attrs.get('Accession')
        if self.is_type:
            if tag == 'GeneList':
                self.is_GeneList = True
            elif tag == 'Gene' and self.ct_gene == 0:
                self.gene_ID = attrs.get('Symbol')
                self.gene_name = attrs.get('FullName')
                self.ct_gene += 1
            elif tag == 'SequenceLocation' and self.is_GeneList == False and self.check_grch == False:
                if attrs.get('Assembly') == 'GRCh38':
                    self.chr = attrs.get('Chr')
                    self.start_num = attrs.get('start')
                    self.stop_num = attrs.get('stop')
                    self.referenceAllele = attrs.get('referenceAlleleVCF')
                    self.alternateAllele = attrs.get('alternateAlleleVCF')
                    self.check_grch = True
            elif tag == 'ProteinExpression' and self.ct_np == 0:
                self.np_acc = attrs.get('sequenceAccessionVersion')
                self.change = attrs.get('change') 
                if self.np_acc and self.np_acc.startswith('NP'):            
                    self.ct_np += 1
            elif tag == 'MolecularConsequence' and self.ct_mc == 0:
                if attrs.get('Type') and 'missense' in attrs.get('Type').lower():
                    self.is_missense = True
                    self.ct_missense += 1
                self.ct_mc += 1
            elif tag == 'Interpretations':
                self.is_interpretations = True
            elif tag == 'Interpretation':
                self.is_interpretation = True
            elif tag == 'Description':
                self.is_description = True
            elif tag == 'DescriptionHistory':
                self.is_desc_hist = True
            
    def end(self, tag):
        if (tag == 'VariationArchive' and self.is_type) or ((self.is_haplotype or self.is_genotype) and tag == 'SimpleAllele'):
            if len(self.intpn) == 1:
                clinical_significance = self.intpn[0].lower()
                if "uncertain" in clinical_significance:
                    self.is_uncertain = True
                    self.ct_uncertain += 1
                elif "conflicting" in clinical_significance:
                    self.is_conflicting = True
                    self.ct_conflicting += 1
                elif "not provided" in clinical_significance:
                    self.is_not_provided = True
                    self.ct_not_provided += 1
            if (self.gene_ID in self.genes_dict.keys()) and self.is_missense and (self.is_uncertain or self.is_conflicting or self.is_not_provided):
                try:
                    self.change = self.change.split('p.')[1]
                    before = aaMapThreeToOne.get(self.change[0:3])
                    after = aaMapThreeToOne.get(self.change[len(self.change) - 3:len(self.change)])
                except: 
                    before = None
                    after = None
                if before and after:  # check if both have a value in aa dict
                    num = self.change[3:len(self.change) - 3]
                    abbreviated_change = before + num + after
                    fasta = ''
                    self.dictlist['gene_ID'].append(self.gene_ID)
                    self.dictlist['gene_name'].append(self.gene_name)
                    self.dictlist['clinical_significance'].append(clinical_significance)
                    self.dictlist['EC_number'].append(self.genes_dict.get(self.gene_ID))
                    self.dictlist['missense_variation'].append(abbreviated_change)
                    self.dictlist['NP_accession'].append(self.np_acc)
                    self.dictlist['ClinVar_accession'].append(self.clinvar_acc)
                    self.dictlist['Chr'].append(self.chr)
                    self.dictlist['start'].append(self.start_num)
                    self.dictlist['stop'].append(self.stop_num)
                    self.dictlist['referenceAllele'].append(self.referenceAllele)
                    self.dictlist['alternateAllele'].append(self.alternateAllele)
            self.ct_gene = 0             
            self.check_grch = False
            self.is_missense = False
            self.is_uncertain = False
            self.is_conflicting = False
            self.is_not_provided = False
            self.ct_np = 0
            self.ct_mc = 0
            self.intpn = []
        if self.is_type:
            if tag == 'GeneList':
                self.is_GeneList = False
            elif tag == 'Interpretations':
                self.is_interpretations = False
            elif tag == 'Interpretation':
                self.is_interpretaion = False
            elif tag == 'Description':
                self.is_description = False
            elif tag == 'DescriptionHistory':
                self.is_desc_hist = False
        if tag == 'VariationArchive':
            self.is_type = False
            self.is_haplotype = False
            self.is_genotype = False
            self.ct +=1
            if self.ct % 10000 == 0:
                print(f'counter: {self.ct}')
                
    def data(self, data):
        if self.is_interpretations and self.is_interpretation and self.is_description and (not self.is_desc_hist):
            self.intpn.append(data)
            
    def close(self):
        print(f"Variations: {self.ct}")
        print(f"Uncertain Significance: {self.ct_uncertain}")
        print(f"Conflicting Report: {self.ct_conflicting}")
        print(f"Not Provided: {self.ct_not_provided}")
        print(f"Missense: {self.ct_missense}")
        print(f"Mutations in the list: {len(self.dictlist['gene_ID'])}")
        print('debug: the file is closed')
        return self.dictlist


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXML(input_path, output_path, gene_dict):
    print('debug: start parcing')
    parser = etree.XMLParser(target=VariationHandler(gene_dict))
    data = etree.parse(input_path, parser)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index = False, header = True)
    return df


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
        self.dictlist = {'pdb_ID': [], 'BLAST_evalue': [], 'hit_from': [], 'hit_to': []}
        self.is_hit_id = False
        self.is_evalue = False
        self.is_hit_from = False
        self.is_hit_to = False
        self.ct_iter = 0
        self.ct = 0
        
    def start(self, tag, attrs):
        if tag == 'Hit':  # there are a few <Hsp> tags within a <Hit> tag
            self.ct_iter += 1
        elif tag == 'Hit_id':
            self.is_hit_id = True
        elif tag == 'Hsp_evalue':
            self.is_evalue = True
        elif tag == 'Hsp_hit-from':
            self.is_hit_from = True
        elif tag == 'Hsp_hit-to':
            self.is_hit_to = True
            
    def end(self, tag):
        if tag == 'Iteration':
            if self.ct_iter == 0:  # when there is no hit, add None to each list
                self.dictlist['pdb_ID'].append(None)
                self.dictlist['BLAST_evalue'].append(None)
                self.dictlist['hit_from'].append(None)
                self.dictlist['hit_to'].append(None)
            self.ct_iter = 0
            self.ct += 1
            if self.ct % 10000 == 0:
                print(f'coutner: {self.ct}')
        elif tag == 'Hit_id':
            self.is_hit_id = False
        elif tag == 'Hsp_evalue':
            self.is_evalue = False
        elif tag == 'Hsp_hit-from':
            self.is_hit_from = False
        elif tag == 'Hsp_hit-to':
            self.is_hit_to = False
        elif tag == 'Hsp':
            self.ct_iter += 1
        
    def data(self, data):
        if self.ct_iter == 1:
            if self.is_hit_id:
                try:
                    pdb = data.split('pdb|')[1]
                except: 
                    print(f'Cannot split pdb {data}')
                    pdb = None
                self.dictlist['pdb_ID'].append(pdb)
            elif self.is_evalue:
                self.dictlist['BLAST_evalue'].append(data)
            elif self.is_hit_from:
                self.dictlist['hit_from'].append(data)
            elif self.is_hit_to:
                self.dictlist['hit_to'].append(data)
    
    def close(self):
        print('Has completed parsing the blast_result xml')
        print(f'Total Count: {self.ct}')
        print(f'Length pdb: {len(self.dictlist["pdb_ID"])}, evalue: {len(self.dictlist["BLAST_evalue"])}, hit_from : {len(self.dictlist["hit_from"])}, hit_to : {len(self.dictlist["hit_to"])}')
        return self.dictlist


def readBlastXML(input_path):
    print('Start parcing')
    parser = etree.XMLParser(target=BlastHandler())
    data = etree.parse(input_path, parser)
    return data    
