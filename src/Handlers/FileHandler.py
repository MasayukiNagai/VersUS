import os
import sys
import configparser

def make_dir(dir_path: str):
    if not os.path.exists(dir_path):
        try:
            os.makedirs(dir_path)
        except OSError:
            sys.exit('Cannot make an output directory.')
    else:
        pass


def checkpath(filepath):
    if not os.path.exists(filepath):
        errmsg = 'Error: Cannot find "{0}". Exiting.'.format(filepath)
        sys.exit(errmsg)


def parse_config(conffile):
    conf = configparser.ConfigParser()
    conf.optionxform = str
    conf.read(conffile, 'UTF-8')
    general_dict = dict(conf.items('general'))
    params_dict = dict(conf.items('parameters'))
    return general_dict, params_dict


def check_config_general(general_dict):
    checkpath(general_dict['genes'])
    checkpath(general_dict['proteomes'])
    checkpath(general_dict['variations'])
    blast = general_dict['blast']
    if blast != 'None':
        checkpath(blast)
    vep = general_dict['vep']
    if vep != 'None':
        checkpath(vep)


def check_config_params(params_dict):
    try:
        int(params_dict['fasta_window'])
    except Exception as e:
        print(f'fastaw_window needs to be an integer {params_dict.get("fasta_window")}')
        print(e)
        sys.exit()
        # raise ValueError(f'fastaw_window needs to be an integer: {params_dict.get("fasta_window")}')
    try:
        float(params_dict['evalue'])
    except Exception as e:
        print(f'evalue needs to be a float value: {params_dict.get("evalue")}')
        print(e)
        sys.exit()
        # raise ValueError(f'evalue needs to be a float value: {params_dict.get("evalue")}')


def format_header(vus_dict):
    header = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
    formatted_header = [item for item in header if item in vus_dict.keys()]
    return tuple(formatted_header)


def write_to_tsv(vus_dict: dict, header: tuple, outfile_path: str):
    with open(outfile_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for vus_id in vus_dict:
            info = []
            for item in header:
                info.append(str(vus_dict[vus_id][item]))
            f.write('\t'.join(info) + '\n')  


def write_to_csv(vus_dict: dict, header: tuple, outfile_path: str):
    with open(outfile_path, 'w') as f:
        f.write(','.join(header) + '\n')
        for vus_id in vus_dict:
            info = []
            keys = vus_dict[vus_id].keys()
            for item in header:
                if item in keys:
                    info.append(str(vus_dict[vus_id][item]))
            f.write(','.join(info) + '\n')


def read_tsv_to_dict(tsv_path: str):
    vus_dict = {}
    with open(tsv_path, 'r') as f:
        header = f.readline().split('\t')
        header = [item.strip() for item in header]
        print(header)
        for i, line in enumerate(f):
            vus_dict[i] = {}
            ls = line.split('\t')
            if (len(header) != len(ls)):
                print(f'The length is different: {line}')
                continue
            for j, item in enumerate(header):
                vus_dict[i][item] = ls[j].strip()
    return vus_dict 


def read_csv_to_dict(csv_path: str):
    vus_dict = {}
    with open(csv_path, 'r') as f:
        header = f.readline().split(',')
        header = [item.strip() for item in header]
        for i, line in enumerate(f):
            vus_dict[i] = {}
            ls = line.split(',')
            for j, item in enumerate(header):
                vus_dict[i][item] = ls[j].strip()
    return vus_dict