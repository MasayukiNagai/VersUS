import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "mnaffairs.intl@gmail.com"

aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


def getHumanGenes(idLists):
    human_enzymes = set()
    protein_info = Entrez.efetch(db='protein', id=idLists, retmode='xml', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
    tree = ET.parse(protein_info)
    root = tree.getroot()
    count = 0  # debug
    for gene in root.findall('./GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier/[GBQualifier_name="gene"]/GBQualifier_value'):
        human_enzymes.add(gene.text)
        count += 1  # debug
    print(count)
    with open('../data/HumanEnzymes.txt', 'w') as filehandle:
        filehandle.seek(0)  # sets the file's current position at 0
        filehandle.truncate()  # sets the logical size of a file to 0, which means to delte its content
        filehandle.writelines('%s\n' % human_enzyme for human_enzyme in human_enzymes)
    return human_enzymes


def getVUS_Clinvar(idLists, genes):
    data = {'Gene': [], 'VUS_protein': [], 'NP': [], 'PDB': []}
    i = 0
    count = 0  # debug
    geneSet = set()
    while i < len(idLists):
        missense_info = Entrez.efetch(db='clinvar', id=idLists[i:i+10000], retmax=10000, rettype='vcv', is_variationid="true", from_esearch="true",
                                  api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
        tree = ET.parse(missense_info)
        root = tree.getroot()
        for variation in root.findall('./VariationArchive'):
            gene = variation.find('./InterpretedRecord/SimpleAllele/GeneList/Gene')
            name = gene.attrib['Symbol']
            geneSet.add(name)
            count += 1
            if name in genes:
                for mutation in variation.findall('.//ProteinExpression'):
                    np_num = mutation.attrib['sequenceAccessionVersion']
                    change = mutation.attrib['change'].split('p.')[1]
                    before = aatranlation.get(change[0:3])
                    after = aatranlation.get(change[len(change) - 3:len(change)])
                    if before and after:  # check if both have a value in aa dict
                        num = change[3:len(change) - 3]
                        abbreviated_change = before + num + after
                        data['Gene'].append(name)
                        data['VUS_protein'].append(abbreviated_change)
                        data['NP'].append(np_num)
                        data['PDB'].append(np.nan)
        if i + 10000 < len(idLists):
            i += 10000
        elif i == len(idLists) - 1:
            break
        else:
            i += len(idLists) % 10000 - 1
        print(count)  # debug
    print(geneSet)
    print(len(geneSet))                
    df = pd.DataFrame(data)
    df.to_csv('MM_enzyme.csv', index = False)
    return df


def getFASTA(np_num, location, beforeMutation, numOfSequence = 5):
    handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text')
    seq_record = SeqIO.read(handle, 'fasta')
    sequence = seq_record.seq
    if location - 1 < len(sequence) and sequence[location - 1] == beforeMutation:
        if 2 * numOfSequence + 1 <= len(sequence):
            if location - 1 - numOfSequence >= 0:
                if location - 1 + numOfSequence < len(sequence):
                    return sequence[location-1-numOfSequence:location+numOfSequence]
                else:
                    return sequence[len(sequence)-1-2*numOfSequence:len(sequence)]
            else:
                return sequence[0:2*numOfSequence+2]
        else:
            return sequence
    else:
        return False
    


# get every enzyme
handleHumanEnzymes = Entrez.esearch(db="protein", retmax=20000,
                                    term="Homo sapiens[Organism] AND RefSeq[Filter] AND (1*[EC/RN Number] OR 2*[EC/RN Number] OR 3*[EC/RN Number] OR 4*[EC/RN Number] OR 5*[EC/RN Number] OR 6*[EC/RN Number])")
readHumanEnzymes = Entrez.read(handleHumanEnzymes)
enzyme_ids = set(readHumanEnzymes['IdList'])
print(len(enzyme_ids))  # check if enzyme_ids have duplicates 
print("{} enzymes are found".format(readHumanEnzymes["Count"]))
# protein_detail = Entrez.efetch(db='protein', id=readHumanEnzymes['IdList'][0], retmode='xml')
# print(readHumanEnzymes['IdList'][0])
# print(protein_detail.read())


# get human genes
human_genes = getHumanGenes(readHumanEnzymes['IdList'])
print(human_genes)
# human_genes = {'ADH1B', 'PPM1B', 'RPP21', 'MARS1', 'KMT5A', 'PPP1CC', 'ALPP', 'PPM1K', 'CBR3', 'MANEAL', 'PPP2CA', 'HYAL4', 'LIPT2', 'SUCLG2', 'UROC1', 'PLA1A', 'G6PC3', 'EHHADH', 'MMUT', 'MTHFD1', 'PXDN', 'P4HB', 'ACACA', 'INPP5D', 'CDC25A', 'QRSL1', 'AMY1B', 'PPIA', 'XPNPEP1', 'PON2', 'HACD3', 'PLCZ1', 'PTPN23', 'PLD3', 'HSD11B1', 'SDSL', 'PPP2CB', 'NSUN3', 'TGDS', 'ELOVL7', 'LDHC', 'PHOSPHO2', 'PLD2', 'AS3MT', 'PGGHG', 'PDE4C', 'AGO2', 'XRN2', 'NAA40', 'HLCS', 'APIP', 'RPAP2', 'TOP1MT', 'ZDHHC16', 'CYP1A2', 'DIMT1', 'CRELD1', 'HADH', 'AMY2A', 'NCEH1', 'PDE6G', 'DCXR', 'MPG', 'AARS2', 'DPH6', 'ESD', 'PPP5C', 'RTCB', 'RANBP2', 'CERS3', 'MDH1B', 'GANAB', 'PPEF2', 'PNMT', 'DNASE2', 'C1QTNF3', 'PPP2R2B', 'HAPLN3', 'PLB1', 'KL', 'NT5C3B', 'UBE2I', 'DNASE1L2', 'ALAD', 'CHIA', 'HOGA1', 'LYZL2', 'METTL15', 'CHDH', 'CNOT6L', 'HAGHL', 'CYP2J2', 'IARS2', 'EXOG', 'SSH3', 'SUCLA2', 'MAN2B2', 'PDE3A', 'MBOAT1', 'ACP1', 'CTDNEP1', 'SEPTIN1', 'GDE1', 'MAN1C1', 'IDUA', 'ZDHHC11', 'IDH1', 'PPIAL4D', 'BST1', 'TPTE2', 'PLA2G2C', 'GGCX', 'MPI', 'HINT3', 'SCLY', 'DNASE1L3', 'FKBP8', 'IDH3B', 'CDO1', 'GNAT3', 'SND1', 'PPIL3', 'PDIA3', 'CTPS2', 'MLYCD', 'KAT8', 'GPI', 'TAF1', 'DSEL', 'PIN4', 'ALOXE3', 'PPM1E', 'PPM1F', 'RPIA', 'GPD1L', 'PAN2', 'PTPN12', 'PPCDC', 'CEL', 'SPOUT1', 'LIPG', 'POP4', 'MAN2A2', 'CCR4', 'NAA80', 'AGPAT4', 'PGLS', 'METTL14', 'GLYATL2', 'THEM5', 'PPIAL4H', 'ACOT2', 'ALAS1', 'PLA2G4E', 'ZDHHC7', 'FAHD2A', 'CES3', 'ABHD15', 'N6AMT1', 'PPP5D1', 'UXS1', 'PPP3CB', 'GLB1L2', 'MGLL', 'PDE2A', 'RPE65', 'RAD9A', 'AOAH', 'GAD1', 'PEMT', 'AADACL3', 'DUSP15', 'RARS1', 'PPP2R1A', 'PFAS', 'HYI', 'PHOSPHO1', 'LIG4', 'NAXE', 'PGAM5', 'CLC', 'NSUN2', 'ADO', 'CES2', 'PTPN3', 'DCUN1D3', 'ATOX1', 'CHAC1', 'DUSP8', 'PGAM2', 'AHCYL1', 'ALOX15B', 'ACSF3', 'HGD', 'ELP3', 'CBR4', 'LDHA', 'ENPP6', 'ARSL', 'CTPS1', 'PLA2G4F', 'FKBP14', 'URAD', 'RDH12', 'HIBCH', 'MTMR7', 'IDS', 'PLA2G2F', 'ELAC2', 'CDY1', 'SAP130', 'MCAT', 'EPRS1', 'PDE3B', 'PDE6C', 'SSH2', 'KAT5', 'CFH', 'EHMT2', 'FBP2', 'PRDX5', 'HSD17B11', 'GPHN', 'OVCA2', 'ACOT1', 'PPM1J', 'ENO1', 'SI', 'MRM3', 'TPP1', 'PTPRC', 'ABHD10', 'PPM1L', 'ENPP2', 'PLA2G7', 'NAA10', 'AUH', 'FUCA1', 'ABHD17C', 'ADPRHL2', 'ODC1', 'BCDIN3D', 'SEPTIN5', 'PTPN6', 'PPIAL4C', 'MAN1A1', 'PRMT1', 'TYW1B', 'ZDHHC22', 'ST20-MTHFS', 'CTBS', 'NARS2', 'SHMT1', 'ACOT7', 'NAGA', 'DHRS11', 'TDP1', 'PCCA', 'PTPN7', 'OGG1', 'GNMT', 'MTMR2', 'RNPEP', 'HNMT', 'ZDHHC8', 'ACSM5', 'KMT2C', 'GADL1', 'ASL', 'MDP1', 'ACAT2', 'CPT1A', 'BCO1', 'FARS2', 'NNMT', 'ERAP2', 'LPCAT3', 'GLYAT', 'SARS1', 'HCCS', 'NSD2', 'RPEL1', 'PNKP', 'CMTR1', 'ALPL', 'HYAL1', 'SPAM1', 'SDR42E1', 'GPX8', 'AKR1C2', 'PTRH2', 'DUSP21', 'PNLIPRP1', 'AGPAT2', 'PLA2G15', 'FAP', 'TREH', 'HMGCLL1', 'ACOT8', 'NCOA3', 'SSH1', 'ACP6', 'PPP1R1B', 'PTPA', 'FKBPL', 'PISD', 'PLA2G12A', 'DNASE2B', 'NOCT', 'CPPED1', 'PIGV', 'PXDNL', 'IDH3G', 'PLA2G4B', 'ACP5', 'PRDX4', 'ECHDC1', 'TYMS', 'TKT', 'KAT6A', 'METTL11B', 'NT5C', 'TKTL2', 'GC', 'DPEP1', 'SMYD1', 'CYP1B1', 'MIF', 'PRDM9', 'HDC', 'PTPN21', 'PDP1', 'SETD7', 'INPP4B', 'PEPD', 'UBA6', 'CDC25B', 'CHAT', 'FKBP3', 'NADSYN1', 'NEIL1', 'SETMAR', 'DPP7', 'MECOM', 'TARS3', 'DTD2', 'NARS1', 'CBSL', 'DNPEP', 'ERI2', 'DPP3', 'SEC16A', 'NKTR', 'CFTR', 'PLAAT2', 'COQ5', 'INTS11', 'DGAT1', 'YARS2', 'GPAT4', 'EZH1', 'SYNJ2', 'LARS2', 'CARS2', 'RNASE3', 'PLCD4', 'SETD1B', 'AMT', 'CDKN3', 'SMYD3', 'ACSL6', 'PDIA4', 'PPP1R3A', 'DNMT1', 'PTPN18', 'EDEM3', 'HTATIP2', 'PTPN11', 'MTMR6', 'PTGES2', 'FKBP10', 'LIPC', 'PTPRJ', 'PDE12', 'MIP', 'NAA30', 'NOL9', 'ABHD16A', 'ASNS', 'UNG', 'RNASEK', 'DUSP6', 'DCT', 'OSGEPL1', 'GPLD1', 'PIP4P2', 'HSD17B3', 'CARS1', 'PNLIP', 'APEX2', 'SPACA5B', 'NAPEPLD', 'ATF2', 'TREX2', 'UGDH', 'SETD2', 'CAMKMT', 'DUSP4', 'PTEN', 'RDH11', 'PNLIPRP2', 'ACAT1', 'HSD17B8', 'NACA', 'UMPS', 'RDH16', 'MOCS1', 'GALNS', 'SETD1A', 'PNPLA4', 'CRAT', 'CARM1', 'AMACR', 'MGST3', 'ZDHHC20', 'DNASE1', 'G6PC2', 'KATNAL2', 'BRD4', 'ENO3', 'AKR1B10', 'MAN1A2', 'PLCE1', 'GPX4', 'AARS1', 'ABHD3', 'PLAAT1', 'PAICS', 'GCAT', 'GLCE', 'PC', 'MARS2', 'ARSJ', 'TRIM25', 'CA5A', 'PHYKPL', 'PRMT6', 'NT5M', 'DNASE1L1', 'PPTC7', 'ADSL', 'CLOCK', 'ILVBL', 'PMM1', 'CBS', 'MBLAC1', 'SARS2', 'HPD', 'CA6', 'FASN', 'ZDHHC9', 'FKBP9', 'DOT1L', 'ACSF2', 'ENPEP', 'MAN2C1', 'PTPRT', 'PPP3CA', 'AASDH', 'AWAT1', 'LPCAT4', 'HEXB', 'SORD', 'SMPD1', 'IDI2', 'DUSP26', 'CERS5', 'ACACB', 'MANEA', 'PUS10', 'D2HGDH', 'ABHD6', 'PUDP', 'PTPRN2', 'NT5C2', 'AGO3', 'NAT9', 'GSTZ1', 'DUSP9', 'ADPRH', 'PPEF1', 'EYA3', 'PTRHD1', 'MCUR1', 'LIG3', 'L3HYPDH', 'KIF3B', 'CYP2S1', 'AHCYL2', 'ADHFE1', 'PLCD1', 'BHMT', 'ZDHHC14', 'CA8', 'ALDOA', 'PRMT7', 'TOP2A', 'EPHX2', 'TRMT10A', 'PTPN1', 'PDE8A', 'GNPAT', 'PPIC', 'SDR16C5', 'EYA2', 'MOGAT3', 'GNE', 'KAT2B', 'ARSD', 'SDR42E2', 'SESN2', 'GAA', 'PTP4A1', 'DSE', 'PLA2G6', 'DPP9', 'ALDOC', 'SETDB2', 'PNPLA6', 'PNLDC1', 'ACSL5', 'PLA2G4A', 'ACSS3', 'FH', 'DHDH', 'DUSP10', 'LCMT1', 'PPP1CA', 'MTMR3', 'LYZL1', 'PRDM5', 'ENOSF1', 'PTPN5', 'GRHPR', 'AWAT2', 'HEMK1', 'TPTE', 'TRMT44', 'PTPRU', 'GATM', 'FTCD', 'SKP1', 'DUSP5', 'GPAT2', 'KYAT1', 'PTPRF', 'ADH5', 'NT5E', 'LGALS13', 'OLAH', 'PLPPR2', 'ISG20', 'HGSNAT', 'MBOAT7', 'PNKD', 'ACOT6', 'ACSS1', 'ENGASE', 'FKBP1B', 'PLA2G3', 'CDY2B', 'RDH8', 'TRMT5', 'PRDX6', 'FKBP6', 'PSPH', 'HPGD', 'RPP30', 'ACSBG1', 'SIRPA', 'AKR1C3', 'CTDSP2', 'PON3', 'PTPRZ1', 'MUS81', 'CA3', 'DNMT3B', 'SAMHD1', 'MAN2B1', 'SLC27A3', 'IDI1', 'RPP38', 'INPP5F', 'CROT', 'NR1H3', 'PLCG2', 'RPS3', 'GPX6', 'NQO2', 'MDH2', 'RPP14', 'PRDX3', 'EXOSC9', 'PCCB', 'ESCO2', 'PRMT9', 'PPP4C', 'NACAD', 'GABPB2', 'HSD11B2', 'DHRS7B', 'PTP4A2', 'PPCS', 'ASH1L', 'PPP2R3A', 'CA4', 'MTFMT', 'VARS1', 'MTHFS', 'PPIB', 'G6PC', 'OTC', 'ALPG', 'PPIAL4G', 'LYG2', 'ANG', 'PDE6A', 'NSDHL', 'GPD2', 'BHMT2', 'PFKFB2', 'ABHD5', 'NT5C1B-RDH14', 'PTPN9', 'CDYL', 'DAGLB', 'HSD17B10', 'LRTOMT', 'CELF6', 'PPP1R12B', 'CEMIP', 'RDH5', 'SESN1', 'MANBA', 'LTA4H', 'HSD17B12', 'PRR29', 'PPIL4', 'ACSM2A', 'LYZL6', 'CD38', 'ACOT12', 'SETD3', 'CNP', 'GLYATL1B', 'VPS29', 'PTPN20', 'GCLM', 'PPM1M', 'RPP40', 'SPR', 'LCT', 'XPNPEP2', 'HARS2', 'HENMT1', 'CA7', 'PHLPP1', 'PPIAL4E', 'PLPP5', 'HAO1', 'AKR1A1', 'HACD4', 'ACOT11', 'SPAST', 'CA9', 'PGAM1', 'PPP1CB', 'ACMSD', 'ACO2', 'IDH2', 'PFKFB1', 'GNS', 'GCLC', 'ENO4', 'TRMT10B', 'EHMT1', 'PLPP3', 'LDHAL6A', 'HSD17B2', 'ALOX12B', 'ARSB', 'POP1', 'HYAL3', 'PORCN', 'TREX1', 'RARS2', 'LYPLA1', 'HARS1', 'HAAO', 'PRMT5', 'PDE4A', 'KATNA1', 'PTGIS', 'POP5', 'GPX7', 'ENPP3', 'DUPD1', 'DICER1', 'PPP6C', 'PLA2G2D', 'PRDM6', 'ECHS1', 'PAFAH2', 'PTP4A3', 'BDH1', 'FUOM', 'DTD1', 'GTF2B', 'HAO2', 'ALKBH1', 'GAMT', 'NEU2', 'GALE', 'ACE', 'PLAAT3', 'NMT1', 'NMT2', 'MINPP1', 'UBE2NL', 'HYAL2', 'CPT2', 'OCLN', 'PRMT2', 'CARNS1', 'PPIH', 'CTDP1', 'IMPDH2', 'ZDHHC17', 'PTGDS', 'PTPRH', 'PPP2R2A', 'DUSP19', 'GGACT', 'DUSP11', 'FKBP1A', 'ARSK', 'MTMR8', 'LPIN1', 'ASS1', 'INMT', 'ETHE1', 'CPS1', 'SDS', 'NAA50', 'NPL', 'SPO11', 'PLPP4', 'KAT7', 'NEIL2', 'TTLL10', 'ACP3', 'DGLUCY', 'GUSB', 'ACSM4', 'GLYATL3', 'METTL1', 'PLCG1', 'POP7', 'RNMT', 'HSD17B4', 'KAT6B', 'IMPA2', 'GLB1', 'BRCC3', 'ME1', 'PFKFB4', 'PTPRE', 'HTD2', 'METAP2', 'THNSL2', 'AMY1C', 'KMT5B', 'FARSB', 'PPT1', 'METTL7A', 'RNASEH1', 'LCLAT1', 'DUSP7', 'PRDM16', 'ELAC1', 'PLA2G10', 'AGL', 'ACSS2', 'PRMT8', 'METTL16', 'MCM3AP', 'ARSI', 'BCHE', 'AKR1E2', 'DUSP14', 'NAT10', 'OXSM', 'KYAT3', 'SAT2', 'GTF3C4', 'METTL6', 'KMT2B', 'DUSP18', 'ACOT4', 'FGFRL1', 'STAMBP', 'PDIA2', 'ALOX5', 'HACD1', 'RPE', 'INPP5E', 'SHMT2', 'ZDHHC4', 'NOTUM', 'ACSL4', 'PLA2G4C', 'FBP1', 'TKTL1', 'ZDHHC6', 'NAT8L', 'SPTLC1', 'CA1', 'CES5A', 'RNASEH2A', 'UBA1', 'TRDMT1', 'NANP', 'ACAA2', 'ZDHHC11B', 'KMT2A', 'NAT1', 'SPTLC2', 'SDR9C7', 'GALC', 'CBR1', 'BUD23', 'HARBI1', 'ELOVL1', 'LIG1', 'ZDHHC19', 'ABHD14B', 'KITLG', 'SPCS1', 'MDH1', 'ETNPPL', 'UBASH3B', 'ADH7', 'PPIF', 'PPIL1', 'PPIE', 'LCMT2', 'MGMT', 'GPX2', 'GPCPD1', 'LNPEP', 'BPGM', 'PLCH1', 'MCCC1', 'AKR1C4', 'GLA', 'MGAM', 'ACSM1', 'DAGLA', 'SMG6', 'DDC', 'BPNT1', 'PDE5A', 'ENDOV', 'PPP3R1', 'EXOSC7', 'DROSHA', 'TALDO1', 'DHRS4', 'PLPP1', 'LIPT1', 'NEU3', 'INPP5A', 'CYP1A1', 'SPTLC3', 'HDDC3', 'PPP3R2', 'PTPN13', 'PDE6H', 'NAGLU', 'ENO2', 'NTHL1', 'CEMIP2', 'CMTR2', 'FTSJ3', 'ZDHHC3', 'LTC4S', 'HAGH', 'GARS1', 'LPL', 'DUSP28', 'GPX5', 'ADH1C', 'IMPA1', 'CNDP2', 'ACSM3', 'RDH10', 'PPT2', 'WARS1', 'DERA', 'INPPL1', 'PLPPR4', 'PTPRN', 'FEN1', 'GLYATL1', 'LPIN3', 'NEIL3', 'GLUL', 'NPEPPS', 'NT5C1A', 'METTL4', 'ZDHHC23', 'EARS2', 'ZDHHC18', 'EYA4', 'HSD17B7', 'MPO', 'SOAT1', 'ICMT', 'ABHD17A', 'INPP5K', 'AMY2B', 'WWOX', 'HAL', 'CDC25C', 'HACD2', 'PCBD2', 'MTMR4', 'GPAM', 'GBA3', 'ASPG', 'ELOVL6', 'TPP2', 'INPP5J', 'PLCB3', 'ZDHHC13', 'PNPLA8', 'NAGS', 'FUCA2', 'PAFAH1B1', 'PARG', 'CDC14B', 'CHIT1', 'EBP', 'ALOX12', 'HSD17B6', 'PTPN14', 'GPD1', 'PTPDC1', 'PPM1A', 'CPT1B', 'BDH2', 'RPP25', 'AACS', 'CRYL1', 'MCEE', 'SMPD3', 'GART', 'VARS2', 'PTPRM', 'PDE1C', 'SGPL1', 'YARS1', 'PUS7L', 'RDH14', 'ACSL3', 'PTPRQ', 'DNMT3A', 'PMM2', 'FKBP7', 'CA2', 'PAM', 'LAP3', 'PARS2', 'SMPD4', 'CTSC', 'HIBADH', 'TDO2', 'MRI1', 'THEM4', 'PTRH1', 'HAT1', 'GNPNAT1', 'PTPRK', 'FKBP5', 'DARS2', 'PGP', 'GPAT3', 'SOAT2', 'NAGPA', 'SRR', 'TARS1', 'HMGCL', 'ANXA3', 'LDHAL6B', 'SSU72', 'DLAT', 'SGSH', 'PCK1', 'GDPD2', 'HPX', 'ALAS2', 'METTL3', 'CCDC155', 'FKBP2', 'DPH5', 'CPT1C', 'AKR1B1', 'STS', 'ZDHHC1', 'STAMBPL1', 'SGSM2', 'ME2', 'METAP1', 'BCO2', 'MYORG', 'DBT', 'TDG', 'AIFM1', 'TRMT13', 'WARS2', 'ACP7', 'SYNJ1', 'ZDHHC2', 'PPP3CC', 'SLC27A4', 'TTL', 'ABHD12', 'ACSL1', 'PUS1', 'PDE7A', 'CTDSPL2', 'PLCB4', 'PPME1', 'LPCAT1', 'ACAA1', 'PPIAL4F', 'PLD4', 'EZH2', 'RIMKLB', 'NEU1', 'GGCT', 'PGM2', 'DPEP2', 'NSD1', 'RDH13', 'PDE6B', 'CTDSPL', 'PTPRB', 'FARSA', 'CLYBL', 'PHGDH', 'PPP2R1B', 'LDHD', 'PPP1R3D', 'AKR1B15', 'KATNAL1', 'IMPAD1', 'CASD1', 'PGAM4', 'DUSP16', 'HEXD', 'PPM1G', 'EXD2', 'FOLR3', 'PTPRO', 'ARSG', 'PDE4B', 'MOGAT2', 'AGPAT5', 'PPIG', 'AHCY', 'SUV39H1', 'MTHFD1L', 'HMGCR', 'APTX', 'PTPRS', 'OGA', 'PGD', 'PDE11A', 'PLPP2', 'PLA2G5', 'ELOVL3', 'PLCB1', 'DPEP3', 'PPWD1', 'PTPRG', 'AANAT', 'APEX1', 'FKBP15', 'LARS1', 'ATIC', 'DPP8', 'MANSC1', 'UROS', 'CMBL', 'ADSS2', 'GMPS', 'SPACA5', 'ODR4', 'FPGS', 'PRDM7', 'CHAC2', 'CDY2A', 'DDT', 'PTGES', 'UROD', 'ALOX15', 'DUSP1', 'NT5C3A', 'LCAT', 'CA12', 'NT5DC2', 'TARS2', 'ALPI', 'LPCAT2', 'PGM1', 'RGN', 'DUSP22', 'MUTYH', 'TMEM86B', 'DUSP2', 'PIR', 'TPI1', 'NAA11', 'DUSP12', 'TRMT61B', 'PRXL2B', 'ECI2', 'SLC27A2', 'PTPRA', 'SETDB1', 'PTPN2', 'TRMT61A', 'AMY1A', 'CA14', 'NSD3', 'PDIA5', 'TDP2', 'CREBBP', 'PRDM2', 'HSD3B1', 'DNAJC6', 'TOP3A', 'HPSE', 'ILKAP', 'DLST', 'INPP5B', 'DDTL', 'EPM2A', 'TSTA3', 'OSGEP', 'ADH4', 'AIP', 'ACSBG2', 'CNOT6', 'CDC14A', 'TOP2B', 'LGALS16', 'CA5B', 'PLA2G12B', 'PTPN4', 'BAAT', 'ALKBH8', 'PARN', 'ATAT1', 'KDSR', 'ACP2', 'FAHD1', 'PPM1N', 'PTPRR', 'RAD1', 'ADI1', 'MRPL58', 'MAN1B1', 'ADH1A', 'LDHB', 'UBE2M', 'ZFAND1', 'AGPAT1', 'AMD1', 'PLCD3', 'ENPP1', 'DUSP23', 'FAN1', 'ACSM6', 'CTH', 'GPX1', 'EYA1', 'OCRL', 'TPO', 'PNPLA3', 'PHLPP2', 'GBA2', 'PRORP', 'PGM3', 'NAT2', 'PLPPR3', 'GANC', 'EXOSC8', 'ECI1', 'SIAE', 'ZDHHC5', 'AGPAT3', 'PCMT1', 'APLF', 'ABCC4', 'CAD', 'MBOAT2', 'PTPRD', 'PPM1H', 'PDE6D', 'DUSP13', 'PAFAH1B3', 'PIN1', 'CHCHD2', 'PLPP7', 'ACP4', 'RNASE4', 'SUCLG1', 'MYSM1', 'EPX', 'SMPD2', 'MVD', 'EPHX3', 'TGS1', 'ZDHHC21', 'USP50', 'PDE4D', 'CNOT8', 'TPMT', 'COMT', 'TOP3B', 'GSS', 'NAT8', 'MTR', 'ALDOB', 'LYPLAL1', 'PLA2G4D', 'PIP4P1', 'LRAT', 'QARS1', 'ADH6', 'DGAT2', 'CA13', 'SAT1', 'MAN2A1', 'LPIN2', 'ENOPH1', 'TOP1', 'ELOVL5', 'HSPD1', 'ZDHHC15', 'MOGAT1', 'DARS1', 'ABHD2', 'PLCH2', 'ZNF596', 'CARNMT1', 'HADHA', 'PTGES3', 'DPH7', 'DIS3L', 'ARSA', 'DPP4', 'DUSP3', 'IDH3A', 'KARS1', 'NDUFAF7', 'DHRS7C', 'ZDHHC24', 'SLC27A5', 'ACSM2B', 'PDE9A', 'ACO1', 'LPO', 'TRMT2B', 'PPM1D', 'NCOA1', 'GPX3', 'PCBD1', 'AKR7A2', 'EPHX1', 'RTCA', 'RPL19', 'PTPMT1', 'PDE8B', 'NAA20', 'AADAC', 'PXYLP1', 'SMYD2', 'MRPL44', 'IARS1', 'CDY1B', 'PPIAL4A', 'NTMT1', 'HSD3B7', 'SUV39H2', 'ENPP7', 'PNPLA2', 'ASMT', 'MGAM2', 'DHRS3', 'HSD17B1', 'NEU4', 'KMT5C', 'RIMKLA', 'MRM2', 'CES1', 'LIPA', 'NFX1', 'H6PD', 'MTMR1', 'PRDX2', 'PPID', 'CERS4', 'ADSS1', 'L2HGDH', 'PDE7B', 'ELOVL2', 'INPP1', 'PAN3', 'COQ3', 'RENBP', 'CAT', 'ARSF', 'PDE1B', 'PDIA6', 'PAFAH1B2', 'G6PD', 'PIGK', 'RPUSD4', 'KMT2D', 'HSD11B1L', 'TRMT10C', 'NAPRT', 'TMX3', 'ZDHHC12', 'PPIL6', 'AKR1C1', 'DHRS9', 'UBLCP1', 'MYH14', 'CNOT7', 'KAT2A', 'NT5C1B', 'PUS3', 'PLA2G1B', 'LIPE', 'CSKMT', 'KMT2E', 'GDPD1', 'ME3', 'XPNPEP3', 'PPP2R2C', 'MCCC2', 'EP300', 'PDE10A', 'CRELD2', 'PLAAT4', 'GBA', 'MTM1', 'RNPEPL1', 'MBOAT4', 'PNLIPRP3', 'CTDSP1', 'TRMT11', 'RNASEL', 'CWC27', 'TMPPE', 'LSS', 'PCIF1', 'METAP1D', 'PCK2', 'TYW1', 'HEXA', 'IMPDH1', 'ACOD1', 'CERS2', 'TIGAR', 'PLCB2', 'ACHE', 'GALM', 'PLA2G2E', 'RNGTT', 'GMDS', 'FKBP4', 'TRMT1', 'MOGS', 'MPPE1', 'LYPLA2', 'PON1', 'FTSJ1', 'HADHB', 'PTPN22', 'PDE1A', 'SPCS2', 'PLD1', 'GDPD5', 'PTS', 'GABPB1', 'ISYNA1', 'ELOVL4', 'INPP4A', 'HPGDS', 'LIPF', 'TYW3', 'PGAP6', 'IDO1', 'TFB1M', 'PNPLA7', 'PDXP', 'ABHD17B', 'PPP1R1A', 'FKBP11', 'IDO2', 'PRDX1', 'ANPEP', 'PLD6', 'PFKFB3', 'PLAAT5', 'PDP2', 'CSAD', 'LYZ', 'PLA2G2A', 'SETD5', 'NAXD', 'PUS7', 'NAA60', 'GAD2', 'HSD3B2', 'SCP2', 'PPP1R2', 'TBXAS1', 'CNDP1'}
# print(len(human_genes))

# get every missense genes
handleMissense = Entrez.esearch(db='clinvar', retmax=200000,
                                term='( ( ("clinsig has conflicts"[Properties]) OR ("clinsig vus"[Properties]) ) AND ("missense variant"[molecular consequence] OR "SO 0001583"[molecular consequence]))')
recordsMissense = Entrez.read(handleMissense)
print("{} variants are found".format(recordsMissense['Count']))
# print(recordsMissense)
# detailMissesnse = Entrez.efetch(db='clinvar', id=recordsMissense['IdList'][14000:14003], rettype='vcv', is_variationid="true", from_esearch="true")
# print(detailMissesnse.read())

# df_Clinvar = getVUS_Clinvar(recordsMissense['IdList'], human_genes)
# print(df_Clinvar)