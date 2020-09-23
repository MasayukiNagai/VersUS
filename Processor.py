import csv
import os
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
        fasta_dict = {}
        for root, d_names, file_names in os.walk(seq_dir_path):
            for filename in file_names:
                fname = os.path.join(root, filename)
                with open(fname, 'r') as f:
                    print('opened a fasta file')
                    np_num = ''
                    sequence = ''
                    for line in f:
                        if line[0] == '>':
                            if sequence != '':
                                fasta_dict[np_num] = sequence
                                np_num = ''
                                sequence = ''                    
                            i = 1
                            while line[i] != ' ':
                                np_num += line[i]
                                i += 1
                        else:
                            line = line.strip('\n')
                            sequence += line
        print(f'The number of items in the fasta dict: {len(fasta_dict)}')
        return fasta_dict


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
    def crop_seq(self, sequence: str, pos: int, ref: int, seq_range: int) -> str:
        if pos - 1 < len(sequence) and sequence[pos - 1] == ref:
            if pos - 1 - seq_range <= 0:
                proteinSeq = sequence[0:2 * seq_range + 1]
            elif pos - 1 + seq_range > len(sequence) - 1:
                proteinSeq = sequence[0 if len(sequence) - 1 - 2 * seq_range <= 0 else len(sequence) - 1 - 2 * seq_range:len(sequence)]
            else:
                proteinSeq = sequence[pos - 1 - seq_range : pos + seq_range]
            return proteinSeq
        else:
            return '' 
    
    # add fasta sequence to vus dict
    # return vus dict and unfound_seq set
    def add_seq_to_dict(self, vus_dict: dict, fasta_dict: dict, seq_range: int):
        unfound_seq = set()
        seq_list = []
        for vus_id in vus_dict.keys():
            ref = vus_dict[vus_id]['ref']
            pos = vus_dict[vus_id]['pos']
            np_num = vus_dict[vus_id]['NP_accession']
            try:
                seq = fasta_dict[np_num]
                seq_cropped = self.cropFASTA(seq, pos, ref, seq_range)
            except:
                seq = ''
                unfound_seq.add(np_num)
            vus_dict[vus_id]['sequence'] = seq_cropped
        return vus_dict, unfound_seq

    
    # make FASTA format text file from dataframe for blast search
    def make_fasta_for_blast(self, vus_dict: dict, outfile_path: str):
        # info_tup = ('NP_accession', 'gene_id', 'gene_name', 'sequence')
        with open(outfile_path, 'w') as f:
            for vus_id in range(0, len(vus_dict)):
                np_num = vus_dict[vus_id]['NP_accession']
                gene_id = vus_dict[vus_id]['gene_id']
                gene_name = vus_dict[vus_id]['gene_name']
                seq = vus_dict[vus_id]['sequence']
                tup = (np_num, gene_id, gene_name)
                header = '>' + '\t'.join(tup) + '\n'
                fasta = seq + '\n'
                f.write(header)
                f.write(fasta)


    def blast_locally(self, fasta_path: str, outfile_path: str, evalue: float=10.0, window_size: int=3):
        start = datetime.datetime.now()
        cmd1 = '../ncbi/blast/bin/'
        cmd2 = 'blastp' + ' '         + '-query ' + fasta_path + ' '         + '-db ' + cmd1 + 'pdbaa' + ' '         + '-evalue ' + str(evalue) + ' '         + '-outfmt ' + '5' + ' '         + '-out ' + outfile_path
        cmd = cmd1 + cmd2
        b_cmd = os.system(cmd)
        end = datetime.datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        print(f'Running BLAST took {c[0]} minutes {c[1]} seconds')
        print(cmd + ' : ran with exit code %d' %b_cmd)
    
    
    def add_blast_results(self, vus_dict, blast_dict):
        if set(vus_dict.keys()) == set(blast_dict.keys()):
            print(f'VUS_dict and Blast_dict have different keys. vus_dict: {len(vus_dict)}, blast_dict: {len(blast_dict)}')
            return None
        for common_id in range(0, len(vus_dict)):
            vus_dict[common_id]['pdb_ID'] = blast_dict[common_id]['pdb_ID']
            vus_dict[common_id]['BLAST_evalue'] = blast_dict[common_id]['BLAST_evalue']
            vus_dict[common_id]['hit_from'] = blast_dict[common_id]['hit_from']
            vus_dict[common_id]['hit_to'] = blast_dict[common_id]['hit_to']
        return blast_dict


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

    
    def cropVepOutput(self, vep_file_path, outfile_path):
        cmd = "cat " + vep_file_path + " | "\
            + "grep -v '##'" + " | "\
            + "sed 's/^#\(.*\)/\\1/'" +  " >| "\
            + outfile_path
        crop_cmd = os.system(cmd)
        print(cmd + ' : ran with exit code %d' %crop_cmd)