# Workflow

## Python

### SeqHandler
```
genes_dict = {gene_id: {'ec': ec, 'uniprot_id': uniprot_id}
```

### ClinVarHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897'
              'gene_id': 'GUSB'
              'gend_name': 'glucuronidase beta
              'clinincal_significance': 'uncertain significance'
              'missense_variation': 'R611W'
              'NP_accession': 'NP_000172.2'
              'chr': '7'
              'start': '65961022'
              'stop': '65961022'
              'referenceAllele': 'G'
              'alternateAllele': 'A'
             }
           }
```

### SeqHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897'
              'gene_id': 'GUSB'
              'gend_name': 'glucuronidase beta
              'clinincal_significance': 'uncertain significance'
              'missense_variation': 'R611W'
              'NP_accession': 'NP_000172.2'
              'chr': '7'
              'start': '65961022'
              'stop': '65961022'
              'referenceAllele': 'G'
              'alternateAllele': 'A'
             }
           }
```

### BLASTHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897',
              'gene_id': 'GUSB',
              ...
              'pdb_ID':' 1BHG|A',
              'BLAST_evalue': '1.54947e-11',
              'hit_from': '579',
              'hit_to': '603'
             }
           }
```

### VEPHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897',
              'gene_id': 'GUSB',
              ...
              'gnomAD_AF': '6.631e-06'
             }
           }
```

### CADDHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897',
              'gene_id': 'GUSB',
              ...
              'CADD_score': '26.2'
             }
           }
```

### PTMHandler
```
vus_dict = {vus_id:
             {'ClinVar_accession': '897',
              'gene_id': 'GUSB',
              ...
              'PTM': '2'
             }
           }
```
