from Bio import Entrez
Entrez.email = "masayukinagai_2022@depauw.edu"

# handle = Entrez.esearch(db="protein", retmax=20000,
#                         term="1*[EC/RN Number] AND Homo sapiens[Organism] AND RefSeq[Filter]")
# records = Entrez.read(handle)
# print("{} computational journals found".format(records["Count"]))
# print(records)

# for r in records["IdList"]:
#     detail = Entrez.efetch(db="clinvar", id=r, rettype="fasta", retmode="text")
#     print(detail.read())


handleMissense = Entrez.esearch(db='clinvar', retmax=20000,
                                term='(g6pd[gene] AND ( ( ("clinsig vus"[Properties]) OR ("clinsig has conflicts"[Properties]) ) AND ("missense variant"[molecular consequence] OR "SO 0001583"[molecular consequence])))')
recordsMissense = Entrez.read(handleMissense)
#print("{} computational journals found".format(recordsMissense['Count']))
#print(recordsMissense)
detailMissesnse = Entrez.efetch(db='clinvar', id=650072, rettype='variation', retmode='text')
print(detailMissesnse.read())

# nuccore xml MM



handle.close()
aatranlation = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
