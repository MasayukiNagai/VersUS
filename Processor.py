import os
import csv
import gzip
from Bio import Entrez
from Bio import SeqIO
import datetime

class Processor:

    # read csv file of gene names of human enzymes
    # return dict{gene_names: EC number}
    def readHumanGenesEC(self, path):
        genes_dict = {}
        dup = set()
        with open(path, 'r') as f:
            f.readline()
            reader = csv.reader(f)
            for row in reader:
                key = row[0]
                if key in genes_dict.keys():
                    dup.add(key)
                else:
                    genes_dict[key] = row[1]
        # print(f'{len(dup)} genes are duplicated: ', dup)  # prints, if any, genes duplicated in the file
        return genes_dict, dup


    # makes dictionary of fasta sequences and np number 
    # returns dict{np_num: fasta sequence}
    def make_seq_dict(self, seq_dir_path):
        seq_dict = {}
        ct_np = 0
        for root, d_names, file_names in os.walk(seq_dir_path):
            for filename in file_names:
                fname = os.path.join(root, filename)
                with open(fname, 'r') as f:
                    print(f'opened a fasta file: {fname}')
                    np_num = ''
                    for line in f:
                        if '>' in line:
                            ct_np += 1
                            try:
                                np_num = line.split()[0].split('>')[1]
                                seq_dict[np_num] = ''
                            except:
                                print(f'Error: {line}')
                                return
                            if line[0] != '>':
                                print(f'Not start with >: {line} in {fname}')
                        else:
                            seq_dict[np_num] += line.strip()
        print(f'The number of items in the seq dict: {len(seq_dict)}')
        print(f'The number of np found in the files: {ct_np}')
        return seq_dict


    # fetches fasta sequences for varinants whose sequences aren't in the imported files
    # returns dict{np_num: fasta sequence}
    def get_seq(self, ls_np):
        fasta_dict = {}
        for np_num in ls_np:
            handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key='2959e9bc88ce27224b70cada27e5b58c6b09')
            seq_record = SeqIO.read(handle, 'fasta')
            sequence = seq_record.seq
            fasta_dict[np_num] = sequence
        return fasta_dict


    # crops fasta sequence
    # returns the sequnece cropped in a specified range
    def crop_seq(self, seq: str, pos: int, ref: str, seq_range: int) -> str:
        if pos-1 < len(seq) and seq[pos-1].upper() == ref.upper():
            if pos-1 - seq_range <= 0:
                proteinSeq = seq[0:2*seq_range+1]
            elif pos-1 + seq_range > len(seq) - 1:
                proteinSeq = seq[len(seq)-1-2*seq_range if len(seq)-1-2*seq_range >= 0 else 0:len(seq)]
            else:
                proteinSeq = seq[pos-1-seq_range:pos+seq_range]
            return proteinSeq
        else:
            print('Error: crop_seq')
            return None
    
    # add fasta sequence to vus dict
    # return vus dict and unfound_seq set
    def add_seq_to_dict(self, vus_dict: dict, seq_dict: dict, seq_range: int=12):
        unfound_seq = set()
        seq_list = []
        for vus_id in vus_dict.keys():
            mutation = vus_dict[vus_id]['missense_variation']
            ref = mutation[0]
            pos = int(mutation[1:-1])
            np_num = vus_dict[vus_id]['NP_accession']
            try:
                seq = seq_dict[np_num]
                # seq_cropped = ''
                seq_cropped = self.crop_seq(seq, pos, ref, seq_range)
            except:
                seq_cropped = ''
                unfound_seq.add(np_num)
            if(vus_id==1):
                print(mutation, ref, pos, np_num, seq)
            vus_dict[vus_id]['FASTA_window'] = seq_cropped
        return vus_dict, unfound_seq

    
    # make FASTA format text file from dataframe for blast search
    def make_fasta_for_blast(self, vus_dict: dict, outfile_path: str):
        # info_tup = ('NP_accession', 'gene_id', 'gene_name', 'FASTA_window')
        with open(outfile_path, 'w') as f:
            for vus_id in range(0, len(vus_dict)):
                np_acc = vus_dict[vus_id]['NP_accession']
                gene_id = vus_dict[vus_id]['gene_id']
                gene_name = vus_dict[vus_id]['gene_name']
                seq = vus_dict[vus_id]['FASTA_window']
                tup = (np_acc, gene_id, gene_name)
                header = '>' + '\t'.join(tup) + '\n'
                fasta = seq + '\n'
                f.write(header)
                f.write(fasta)


    def blast_locally(self, fasta_path: str, outfile_path: str, evalue: float=10.0, window_size: int=3):
        start = datetime.datetime.now()
        cmd1 = './ncbi/blast/bin/'  # be careful for this path
        cmd2 = 'blastp' + ' ' \
             + '-query ' + fasta_path + ' '\
             + '-db ' + cmd1 + 'pdbaa' + ' '\
             + '-evalue ' + str(evalue) + ' '\
             + '-outfmt ' + '5' + ' '\
             + '-out ' + outfile_path
        cmd = cmd1 + cmd2
        b_cmd = os.system(cmd)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        print(f'Running BLAST took {c[0]} minutes {c[1]} seconds')
        print(cmd + ' : ran with exit code %d' %b_cmd)
    
    
    def add_blast_results(self, vus_dict, blast_dict):
        if vus_dict.keys() != blast_dict.keys():
            print(f'VUS_dict and Blast_dict have different keys. vus_dict: {len(vus_dict)}, blast_dict: {len(blast_dict)}')
            return None
        for common_id in range(0, len(vus_dict)):
            vus_dict[common_id]['pdb_ID'] = blast_dict[common_id]['pdb_ID']
            vus_dict[common_id]['BLAST_evalue'] = blast_dict[common_id]['BLAST_evalue']
            vus_dict[common_id]['hit_from'] = blast_dict[common_id]['hit_from']
            vus_dict[common_id]['hit_to'] = blast_dict[common_id]['hit_to']
        return vus_dict


    def make_tsv_for_vep(self, vus_dict, outfile_path):
        with open(outfile_path, 'w') as f:
            for vus_id in range(0, len(vus_dict)):
                info = (vus_dict[vus_id]['Chr'], vus_dict[vus_id]['start'], vus_dict[vus_id]['stop'], vus_dict[vus_id]['referenceAllele'], vus_dict[vus_id]['alternateAllele'])
                f.write('\t'.join(info) + '\n')
        # df_sub = df[['Chr', 'start', 'stop', 'referenceAllele', 'alternateAllele']]
        # df_sub.to_csv(output_path, sep='\t', index=False, header=False)
        # return df_sub
    

    def vep_locally(self, vep_input_path, outfile_path):
        start = datetime.datetime.now()  # for counting time necessary to run vep
        cmd1 = '../../ensembl-vep/'  # specify directory
        cmd2 = './vep' + ' '\
             + '-i ' + vep_input_path + ' '\
             + '-o ' + outfile_path + ' '\
             + '--cache ' + '--af_gnomad '\
             + '--appris ' + '--canonical '\
             + '--sift p ' + '--polyphen p '\
             + '--no_check_variants_order '\
             + '--flag_pick ' + '--tab '\
             + '--force_overwrite'
        cmd = cmd1 + cmd2
        vep_cmd = os.system(cmd)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        print(cmd + ' : ran with exit code %d' %vep_cmd)
        print(f'Running Vep took {c[0]} minutes {c[1]} seconds')

    
    def crop_vep_output(self, vep_file_path, outfile_path):
        cmd = "cat " + vep_file_path + " | "\
            + "grep -v '##'" + " | "\
            + "sed 's/^#\(.*\)/\\1/'" +  " >| "\
            + outfile_path
        crop_cmd = os.system(cmd)
        print(cmd + ' : ran with exit code %d' %crop_cmd)

    
    def make_tsv_for_CADD(self, vus_dict, outfile_path):
        header_info = ('CHROM', 'POS', 'ID', 'REF', 'ALT')
        header = '#' + '\t'.join(header_info)
        with open(outfile_path, 'w') as f:
            f.write(header + '\n')
            for vus_id in range(len(vus_dict)):
                chrom = vus_dict[vus_id]['chr']
                pos = str(vus_dict[vus_id]['start'])
                ref = vus_dict[vus_id]['referenceAllele']
                alt = vus_dict[vus_id]['alternateAllele']
                info_tup = (chrom, pos, str(vus_id), ref, alt)
                info = '\t'.join(info_tup)
                f.write(info + '\n')
    

    def read_CADD_results(self, caddfile_path: str) -> dict:
        cadd_dict = {}
        with gzip.open(caddfile_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                ls = line.split('\t')
                chrom = ls[0].strip()
                pos = ls[1].strip()
                ref = ls[2].strip()
                alt = ls[3].strip()
                phred = ls[5].strip()
                key = ref + pos + alt
                if chrom not in cadd_dict.keys():
                    cadd_dict[chrom] = {}
                cadd_dict[chrom][key] = phred
        length = 0
        for chrom in cadd_dict.keys():
            length += len(cadd_dict[chrom].keys())
        print(f'Length of CADD dict: {length}')
        return cadd_dict

    def add_cadd_results(self, vus_dict: dict, cadd_dict: dict):
        unfound_cadd = {}
        for vus_id in vus_dict:
            chrom = vus_dict[vus_id]['chr']
            ref = vus_dict[vus_id]['referenceAllele']
            pos = vus_dict[vus_id]['start']
            alt = vus_dict[vus_id]['alternateAllele']
            key = ref + pos + alt
            try:
                cadd_score = cadd_dict[chrom][key]
            except:
                c_acc = vus_dict[vus_id]['ClinVar_accession']
                unfound_cadd[c_acc] = key
                cadd_score = None
            vus_dict[vus_id]['CADD_score'] = cadd_score
        return vus_dict, unfound_cadd
            

    # def write_to_csv(self, vus_dict: dict, header: tuple, outfile_path: str):
    #     with open(outfile_path, 'w') as f:
    #         f.write(','.join(header) + '\n')
    #         for vus_id in vus_dict:
    #             info = []
    #             for item in header:
    #                 info.append(str(vus_dict[vus_id][item]))
    #             f.write(','.join(info) + '\n')


    def write_to_tsv(self, vus_dict: dict, header: tuple, outfile_path: str):
        with open(outfile_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for vus_id in vus_dict:
                info = []
                for item in header:
                    info.append(str(vus_dict[vus_id][item]))
                f.write('\t'.join(info) + '\n')


    # def read_csv_to_dict(self, csv_path: str):
    #     vus_dict = {}
    #     with open(csv_path, 'r') as f:
    #         header = f.readline().split(',')
    #         header = [item.strip() for item in header]
    #         for i, line in enumerate(f):
    #             vus_dict[i] = {}
    #             ls = line.split(',')
    #             for j, item in enumerate(header):
    #                 vus_dict[i][item] = ls[j].strip()
    #     return vus_dict


    def read_tsv_to_dict(self, tsv_path: str):
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
