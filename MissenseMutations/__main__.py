import pandas as pd
import xml.etree.ElementTree as ET
from Bio import Entrez


Entrez.email = "masayukinagai_2022@depauw.edu"

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}

handleHumanEnzymes = Entrez.esearch(db="protein", retmax=20000,
                                    term="Homo sapiens[Organism] AND RefSeq[Filter] AND (1*[EC/RN Number] OR 2*[EC/RN Number] OR 3*[EC/RN Number] OR 4*[EC/RN Number] OR 5*[EC/RN Number] OR 6*[EC/RN Number])")
readHumanEnzymes = Entrez.read(handleHumanEnzymes)
print("{} enzymes are found".format(readHumanEnzymes["Count"]))
protein_detail = Entrez.efetch(db='protein', id=readHumanEnzymes['IdList'][0], retmode='xml')
# print(readHumanEnzymes['IdList'][0])
# print(protein_detail.read())


def getHumanGenes(idLists):
    human_enzymes = set()

    for i in idLists:
        protein_info = Entrez.efetch(db='protein', id=i, retmode='xml')
        tree = ET.parse(protein_info)
        root = tree.getroot()
        gene = root.find('./GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier/[GBQualifier_name="gene"]/GBQualifier_value')
        human_enzymes.add(gene.text)

    return human_enzymes


human_genes = getHumanGenes(readHumanEnzymes['IdList'][0:150])
print(human_genes)

handleMissense = Entrez.esearch(db='clinvar', retmax=200000,
                                term='( ( ("clinsig has conflicts"[Properties]) OR ("clinsig vus"[Properties]) ) AND ("missense variant"[molecular consequence] OR "SO 0001583"[molecular consequence]))')
recordsMissense = Entrez.read(handleMissense)
print("{} genes are found".format(recordsMissense['Count']))
# print(recordsMissense)
detailMissesnse = Entrez.efetch(db='clinvar', id=558834, rettype='vcv', is_variationid="true", from_esearch="true")
# print(detailMissesnse.read())

def getVUS_Clinvar(idLists, genes):
    data = {'Gene': [], 'VUS_protein': [], 'NP': [], 'PDB': []}
    for i in idLists:
        missense_info = Entrez.efetch(db='clinvar', id=i, rettype='vcv', is_variationid="true", from_esearch="true")
        tree = ET.parse(missense_info)
        root = tree.getroot()
        name = root.find('./VariationArchive/InterpretedRecord/SimpleAllele/GeneList/Gene').attrib['Symbol']
        if name in genes:
            for mutation in root.iter('ProteinExpression'):
                data['Gene'].append(name)
                np = mutation.attrib['sequenceAccessionVersion']
                change = mutation.attrib['change'].split('p.')[1]
                before = aatranlation[change[0:3]]
                after = aatranlation[change[len(change) - 3:len(change)]]
                num = change[3:len(change) - 3]
                abbreviated_change = before + num + after
                data['VUS_protein'].append(abbreviated_change)
                data['NP'].append(np)
    df = pd.DataFrame(data)
    return df


df_Clinvar = getVUS_Clinvar(recordsMissense['IdList'][0:100], human_genes)
print(df_Clinvar)

# def getGenes(xmlfile):
#     tree = ET.parse(xmlfile)
#     root = tree.getroot()
#     qualifierNames = dict()
#     for child in root.iter('GBQualifier'):
#         name = child.find('GBQualifier_name').text
#         qualifierNames[name] = True
#     if 'gene' in qualifierNames:
#         return True
#     else:
#         return False
#
#
# def getHumanEnzymes(idLists):
#     human_enzymes = set()
#
#     for i in idLists:
#         protein_info = Entrez.efetch(db='protein', id=i, retmode='xml')
#         gene = getGenes(protein_info)
#         print(gene)  # delete later
#         if gene:
#             human_enzymes.add(i)
#         protein_info.close()
#         print(i)  # delete later
#
#     return human_enzymes


# human_enzymes = getHumanEnzymes(readHumanEnzymes['IdList'][100:110])
# print(human_enzymes)