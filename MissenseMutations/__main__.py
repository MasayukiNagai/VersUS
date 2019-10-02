import xml.etree.ElementTree as ET
from Bio import Entrez


Entrez.email = "masayukinagai_2022@depauw.edu"

handleHumanEnzymes = Entrez.esearch(db="protein", retmax=20000,
                                    term="Homo sapiens[Organism] AND RefSeq[Filter] AND (1*[EC/RN Number] OR 2*[EC/RN Number] OR 3*[EC/RN Number] OR 4*[EC/RN Number] OR 5*[EC/RN Number] OR 6*[EC/RN Number])")
readHumanEnzymes = Entrez.read(handleHumanEnzymes)
print("{} enzymes are found".format(readHumanEnzymes["Count"]))


def getGenes(xmlfile):
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    qualifierNames = dict()
    for child in root.iter('GBQualifier'):
        name = child.find('GBQualifier_name').text
        qualifierNames[name] = True
    if 'gene' in qualifierNames:
        return True
    else:
        return False


def getHumanEnzymes(idLists):
    human_enzymes = set()

    for i in idLists:
        protein_info = Entrez.efetch(db='protein', id=i, retmode='xml')
        gene = getGenes(protein_info)
        print(gene)  # delete later
        if gene:
            human_enzymes.add(i)
        protein_info.close()
        print(i)  # delete later

    return human_enzymes


# human_enzymes = getHumanEnzymes(readHumanEnzymes['IdList'][100:110])
# print(human_enzymes)

handleMissense = Entrez.esearch(db='clinvar', retmax=20000,
                                term='( ( ("clinsig has conflicts"[Properties]) OR ("clinsig vus"[Properties]) ) AND ("missense variant"[molecular consequence] OR "SO 0001583"[molecular consequence]))')
recordsMissense = Entrez.read(handleMissense)
print("{} genes are found".format(recordsMissense['Count']))
print(recordsMissense)
detailMissesnse = Entrez.efetch(db='clinvar', id=558834, rettype='vcv', is_variationid="true", from_esearch="true")
print(detailMissesnse.read())

# nuccore xml MM
aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
