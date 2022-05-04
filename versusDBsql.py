import mysql.connector
import argparse
import configparser
import os
import sys
import gzip
from Bio import SeqIO, Entrez

aaMapOneToThree = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                   'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                   'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                   'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}

cinvar_link_variation = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/'

class DataBaseEditor:
    def __init__(self):
        self.argument_parser()
        self.gene_table = 'Gene'
        self.mutation_table = 'Mutation'
        self.ec_table = 'Enzyme_class'
        self.fasta_table = 'Fasta'
        # self.clinical_significance = 'Clinical_Significane'
        # self.pdb = 'PDB'
        self.gene_id = 0
        self.mutation_id = 0

    def argument_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--config', '-c', type=str, required=True, dest='conf',
            help='Required; Specify a config file (txt).')
        parser.add_argument(
            '--input', '-i', type=str, required=True, dest='vus',
            help='Required; Specify a VUS file (tsv).')
        parser.add_argument(
            '--ec', '-e', type=str, required=True, dest='ec',
            help='Required; Specify an ec file (txt).')
        parser.add_argument(
            '--fasta', '-f', type=str, required=True, dest='fasta',
            help='Required; Specify a direcotry path to fasta.')
        self.args = parser.parse_args()
        self.vus = self.args.vus
        self.fasta = self.args.fasta
        self.ec = self.args.ec

    def check_infiles(self):
        for f in [self.vus, self.fasta, self.ec]:
            if not os.path.exists(f):
                sys.exit(f'Error: Cannot find "{f}"')

    def parse_config(self):
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(self.args.conf, 'UTF-8')
        self.settings = dict(conf.items('settings'))
        self.host = self.settings['host']
        self.user = self.settings['user']
        self.passwd = self.settings['passwd']
        self.dbname = self.settings['dbname']
        self.apikey = self.settings['apikey']
        Entrez.email = self.settings['email']

    def prepare(self):
        self.cnx = mysql.connector.connect(
            host=self.host,
            user=self.user,
            password=self.passwd,
            database=self.dbname
        )
        self.cur = self.cnx.cursor()

    def close(self):
        self.cur.close()
        self.cnx.close()

    def create_table(self, table_name, name_type_dict):
        query1 = 'CREATE TABLE IF NOT EXISTS ' + table_name + ' ('
        bodylist = []
        for name in name_type_dict.keys():
            bodylist.append(name + ' ' + name_type_dict[name])
        body = ','.join(bodylist)
        query2 = ');'
        query = query1 + body + query2
        self.cur.execute(query)

    def create_gene_table(self):
        nn = ' NOT NULL'
        name_type_dict = {'gene_id': 'SERIAL PRIMARY KEY',
                          'gene_symbol': 'VARCHAR(10)' + nn,
                          'gene_full_name': 'VARCHAR(255)' + nn,
                          'uniprot_id': 'VARCHAR(20)' + nn,
                          'chrom': 'VARCHAR(20)' + nn,
                          'EC_number': 'VARCHAR(100)' + nn,
                          'ec_1': 'INT' + nn,
                          'ec_2': 'INT', 'ec_3': 'INT', 'ec_4': 'INT',
                          'EC_number2': 'VARCHAR(100)',
                          'ec2_1': 'INT',
                          'ec2_2': 'INT', 'ec2_3': 'INT', 'ec2_4': 'INT'}
        self.create_table(self.gene_table, name_type_dict)

    def create_mutation_table(self):
        nn = ' NOT NULL'
        name_type_dict = {'mutation_id': 'SERIAL PRIMARY KEY',
                          'gene_id': 'INT' + nn,
                          'ref': 'CHAR(3)' + nn,
                          'pos': 'INT' + nn,
                          'alt': 'CHAR(3)' + nn,
                          'accession': 'VARCHAR(255)' + nn,
                          'clinical_significance': 'VARCHAR(255)' + nn,
                          'CADD_score': 'FLOAT',
                          'gnomAD_AF': 'FLOAT',
                          'pdb': 'VARCHAR(20)',
                          'fasta_id': 'INT'}  # nn?
        self.create_table(self.mutation_table, name_type_dict)

    def create_ec_table(self):
        nn = ' NOT NULL'
        name_type_dict = {'ec_id': 'SERIAL PRIMARY KEY',
                          'ec_number': 'VARCHAR(20)' + nn,
                          'description': 'VARCHAR(500)' + nn,
                          'class': 'INT ' + nn,
                          'ec_1': 'INT ' + nn,
                          'ec_2': 'INT', 'ec_3': 'INT', 'ec_4': 'INT'}
        self.create_table(self.ec_table, name_type_dict)

    def create_fasta_table(self):
        nn = ' NOT NULL'
        name_type_dict = {'fasta_id': 'SERIAL PRIMARY KEY',
                          'NP_accession': 'VARCHAR(20)' + nn,
                          'fasta': 'TEXT' + nn}
        self.create_table(self.fasta_table, name_type_dict)

    def create_all_tables(self):
        # connect to mysql database
        self.prepare()
        # create tables
        self.create_gene_table()
        self.create_mutation_table()
        self.create_ec_table()
        self.create_fasta_table()
        # close connection
        self.close()

    def drop_table(self, table_name):
        query = 'DROP TABLE IF EXISTS ' + table_name
        self.cur.execute(query)

    def drop_all_tables(self):
        # connect to mysql database
        self.prepare()
        # create tables
        self.drop_table(self.gene_table)
        self.drop_table(self.mutation_table)
        self.drop_table(self.ec_table)
        self.drop_table(self.fasta_table)
        print('---------- Dropped all tables ----------')
        # close connection
        self.close()

    def insert_items(self, table, keytup, values_tuplist):
        query1 = 'INSERT INTO ' + table + ' '
        body1 = '(' + ', '.join(keytup) + ') VALUES '
        # values_query_list = []
        # for values in values_tuplist:
        #     val_q = '(' + ', '.join(values) + ')'
        #     values_query_list.append(val_q)
        params = ['%s'] * len(keytup)
        body2 = '(' + ', '.join(params) + ')'
        query = query1 + body1 + body2
        print(query)
        # print(values_tuplist)
        self.cur.executemany(query, values_tuplist)
        self.cnx.commit()

    def str_editor(self, strings):
        # if strings == 'NULL' or strings == 'null':
        #     return strings
        # if type(strings) is str:
        #     return '"' + strings + '"'
        if type(strings) is str:
            return strings
        elif strings is None:
            return None
        else:
            return str(strings)

    def get_tuplist_for_gene(self, keytup, gene_dict):
        tuplist = []
        for gene_id in gene_dict.keys():
            tup = tuple([self.str_editor(gene_dict[gene_id][key]) for key in keytup])
            tuplist.append(tup)
        return tuplist

    def get_tuplist_for_mut(self, keytup, mutation_dict):
        tuplist = []
        for gene_id in mutation_dict.keys():
            for mutation in mutation_dict[gene_id]:
                tup = tuple([self.str_editor(mutation[key]) for key in keytup])
                tuplist.append(tup)
        return tuplist

    def get_tuplist_for_ec(self, keytup, ec_dict):
        tuplist = []
        for ec_number, item in ec_dict.items():
            tup = tuple([ec_number] + [self.str_editor(item[key]) for key in keytup[1:]])
            tuplist.append(tup)
        return tuplist

    def get_tuplist_for_fasta(self, keytup, fasta_dict):
        tuplist = []
        for np_accession, fasta in fasta_dict.items():
            tup = tuple([np_accession, self.str_editor(fasta)])
            tuplist.append(tup)
        return tuplist

    def register_entries(self):
        self.prepare()
        self.register_fasta(self.fasta, self.vus)
        self.register_vus(self.vus)
        self.register_ec(self.ec)
        self.close()

    def register_vus(self, vus_tsv):
        vus_dictlist = self.parse_vus_data(vus_tsv)
        self.register_genes(vus_dictlist)
        self.add_index_gene()
        self.register_mutations(vus_dictlist)
        self.add_index_mutation()

    def register_genes(self, vus_dictlist):
        gene_dict = dict()
        for vus_dict in vus_dictlist:
            if vus_dict['gene_id'] not in gene_dict.keys():
                gene_dict[vus_dict['gene_id']] = {'gene_symbol': vus_dict['gene_id'],
                                                  'gene_full_name': vus_dict['gene_name'],
                                                  'uniprot_id': vus_dict['uniprot_id'],
                                                  'chrom': vus_dict['chr']}
                ec_dict = self.get_ec_info(vus_dict['EC_number'])
                gene_dict[vus_dict['gene_id']].update(ec_dict)
        keytup = ('gene_symbol', 'gene_full_name', 'uniprot_id', 'chrom',
                  'EC_number', 'ec_1', 'ec_2', 'ec_3', 'ec_4',
                  'EC_number2', 'ec2_1', 'ec2_2', 'ec2_3', 'ec2_4')
        gene_values = self.get_tuplist_for_gene(keytup, gene_dict)
        self.insert_items(self.gene_table, keytup, gene_values)

    def get_ec_info(self, ec_number):
        ec_numbers = ec_number.split(';')
        ec1 = ec_numbers[0]
        ec1_levels = ec1.split('.')
        assert len(ec1_levels) == 4, ec1_levels
        ec_1 = int(ec1_levels[0])
        ec_2 = int(ec1_levels[1]) if self.is_integer(ec1_levels[1]) else -1
        ec_3 = int(ec1_levels[2]) if self.is_integer(ec1_levels[2]) else -1
        ec_4 = int(ec1_levels[3]) if self.is_integer(ec1_levels[3]) else -1
        if len(ec_numbers) > 1:
            ec2 = ec_numbers[1].strip()
            ec2_levels = ec2.split('.')
            assert len(ec2_levels) == 4, ec2_levels
            ec2_1 = int(ec2_levels[0])
            ec2_2 = int(ec2_levels[1]) if self.is_integer(ec2_levels[1]) else -1
            ec2_3 = int(ec2_levels[2]) if self.is_integer(ec2_levels[2]) else -1
            ec2_4 = int(ec2_levels[3]) if self.is_integer(ec2_levels[3]) else -1
        else:
            ec2, ec2_1, ec2_2, ec2_3, ec2_4 = None, None, None, None, None
        ec_dict = {'EC_number': ec1, 'ec_1': ec_1, 'ec_2': ec_2, 'ec_3': ec_3, 'ec_4': ec_4,
                   'EC_number2': ec2, 'ec2_1': ec2_1, 'ec2_2': ec2_2, 'ec2_3': ec2_3, 'ec2_4': ec2_4}
        return ec_dict

    def is_integer(self, n):
        try:
            int(n)
        except ValueError:
            return False
        else:
            return True

    def register_mutations(self, vus_dictlist):
        mutation_dict = dict()
        for vus_dict in vus_dictlist:
            if vus_dict['gene_id'] not in mutation_dict.keys():
                mutation_dict[vus_dict['gene_id']] = []
            mut_entry = self.make_mut_dict(vus_dict)
            mutation_dict[vus_dict['gene_id']].append(mut_entry)
        keytup = ('gene_id', 'ref', 'pos', 'alt', 'accession',
                  'clinical_significance', 'CADD_score', 'gnomAD_AF', 'pdb',
                  'fasta_id')
        vus_values = self.get_tuplist_for_mut(keytup, mutation_dict)
        self.insert_items(self.mutation_table, keytup, vus_values)

    def make_mut_dict(self, vus_dict):
        gene_id = self.get_gene_id(vus_dict['gene_id'])
        fasta_id = self.get_fasta_id(vus_dict['NP_accession'])
        ref_pos_alt = vus_dict['missense_variation']
        # ref_pos_alt = aaMapOneToThree(vus_dict['ref']) + vus_dict['pos'] +  aaMapOneToThree(vus_dict['alt'])
        ref = ref_pos_alt[:3]
        pos = int(ref_pos_alt[3:-3])
        alt = ref_pos_alt[-3:]
        accession = vus_dict['ClinVar_accession']
        mut = {'gene_id': gene_id,
               'ref': ref,
               'pos': pos,
               'alt': alt,
               'accession': accession,
               'clinical_significance': vus_dict['clinical_significance'],
               'CADD_score': vus_dict['CADD_score'],
               'gnomAD_AF': vus_dict['gnomAD_AF'],
               'pdb': vus_dict['pdb_ID'],
               'fasta_id': fasta_id}
        return mut

    def get_gene_id(self, gene_name):
        query = 'SELECT gene_id FROM ' + self.gene_table\
               + ' WHERE gene_symbol = ' + '%s'
        self.cur.execute(query, (gene_name,))
        gene_id_list = self.cur.fetchall()
        assert len(gene_id_list) == 1, \
            f'"{gene_name}" does not exist in Gene table or the Table may have duplicates'
        gene_id = gene_id_list[0][0]
        return gene_id

    def get_fasta_id(self, NP_accession):
        query = 'SELECT fasta_id FROM ' + self.fasta_table\
               + ' WHERE NP_accession = ' + '%s'
        self.cur.execute(query, (NP_accession,))
        fasta_ids = self.cur.fetchall()
        assert len(fasta_ids) <= 1,\
            f'"{NP_accession}" may have duplicates.'
        if len(fasta_ids) == 1:
            fasta_id = fasta_ids[0][0]
        else:
            print(f'"{NP_accession}" does not have a corresponding seq.')
            fasta_id = None
        return fasta_id

    def register_ec(self, ec_file):
        ec_dict = self.parse_EC_numbers(ec_file)
        self.register_ec_numbers(ec_dict)

    def register_ec_numbers(self, ec_dict):
        keytup = ('ec_number', 'description', 'class', 'ec_1', 'ec_2', 'ec_3', 'ec_4')
        ec_values = self.get_tuplist_for_ec(keytup, ec_dict)
        self.insert_items(self.ec_table, keytup, ec_values)

    def register_fasta(self, fasta_dirpath, vus_tsv):
        fasta_dict = self.register_fasta_from_files(fasta_dirpath)
        self.register_fasta_from_fetched(fasta_dict, vus_tsv)
        self.add_index_fasta()

    def register_fasta_from_files(self, fasta_dirpath):
        fasta_dict = self.parse_fasta(fasta_dirpath)
        keytup = ('NP_accession', 'fasta')
        fasta_values = self.get_tuplist_for_fasta(keytup, fasta_dict)
        n = 10000
        chunks=[fasta_values[i:i + n] for i in range(0, len(fasta_values), n)]
        for chunk in chunks:
            try:
                self.insert_items(self.fasta_table, keytup, chunk)
            except:
                print('Exception occured')
                self.close()
                self.prepare()
                self.insert_items(self.fasta_table, keytup, chunk)
        return fasta_dict

    def register_fasta_from_fetched(self, fasta_dict, vus_tsv):
        np_set = set()
        with open(vus_tsv, 'r') as f:
            header = f.readline().rstrip().split('\t')
            # keytup = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'uniprot_id', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
            # assert len(keys) == len(keytup)
            assert header[6] == 'NP_accession',\
                'Error: Confirm the column for "NP_accession"'
            for line in f:
                data = line.rstrip().split('\t')
                np_acc = data[6]
                np_set.add(np_acc)
        unfound_nps = np_set - set(fasta_dict.keys())
        print(f'"{len(unfound_nps)}" unfound NP_accessions: {unfound_nps}')
        comp_fasta_dict = self.fetch_seq(unfound_nps)
        keytup = ('NP_accession', 'fasta')
        fasta_values = self.get_tuplist_for_fasta(keytup, comp_fasta_dict)
        self.insert_items(self.fasta_table, keytup, fasta_values)

    def add_index_gene(self):
        query = f'ALTER TABLE {self.gene_table} ADD INDEX name_id_idx(gene_symbol, gene_id);'
        self.cur.execute(query)

    def add_index_mutation(self):
        query = f'ALTER TABLE {self.mutation_table} ADD INDEX gene_cadd_idx(gene_id, CADD_score DESC);'
        self.cur.execute(query)

    def add_index_fasta(self):
        query = f'ALTER TABLE {self.fasta_table} ADD INDEX np_fasta_idx(NP_accession, fasta_id);'
        self.cur.execute(query)

    def parse_vus_data(self, vus_tsv):
        with open(vus_tsv, 'r') as f:
            header = f.readline().rstrip()
            keys = header.split('\t')
            keytup = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'uniprot_id', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
            assert len(keys) == len(keytup)
            vus_dictlist = []
            for line in f:
                data = line.rstrip().split('\t')
                assert len(data) == len(keytup)
                vus_dict = dict(zip(keytup, data))
                for key in keytup:
                    if vus_dict[key] == 'None':
                        vus_dict[key] = None
                    elif key == 'missense_variation':
                        vus_dict[key] = self.map_aa_one_to_three(vus_dict[key])
                vus_dictlist.append(vus_dict)
        return vus_dictlist

    def map_aa_one_to_three(self, missense):
        ref = aaMapOneToThree[missense[0]]
        alt = aaMapOneToThree[missense[-1]]
        pos = missense[1:-1]
        return ref + pos + alt

    def parse_EC_numbers(self, ec_file):
        ec_dict = {}
        with open(ec_file, 'r') as f:
            for line in f:
                data = line.rstrip().split()
                ec = data[0]
                desc = ' '.join(data[1:])
                if "undefined" in desc:
                    continue
                ec_levels = ec.split('.')
                class_level = len(ec_levels)
                ec_1 = int(ec_levels[0])
                ec_2 = int(ec_levels[1]) if class_level > 1 else None
                ec_3 = int(ec_levels[2]) if class_level > 2 else None
                ec_4 = int(ec_levels[3]) if class_level > 3 else None
                ec_dict[ec] = {'description': desc, 'class': class_level, 'ec_1': ec_1, 'ec_2': ec_2, 'ec_3': ec_3, 'ec_4': ec_4}
        return ec_dict

    def parse_fasta(self, fasta_dirpath):
        fasta_dict = {}
        for root, _, file_names in os.walk(fasta_dirpath):
            for filename in file_names:
                fname = os.path.join(root, filename)
                with gzip.open(fname, 'rt') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        fasta_dict[record.id] = str(record.seq)
        print(f'Finish processing {len(fasta_dict)} sequences.')
        return fasta_dict

    def fetch_seq(self, np_ls):
        fasta_dict = {}
        for np_num in np_ls:
            handle = Entrez.efetch(db='protein', id=np_num, rettype='fasta', retmode='text', api_key=self.apikey)
            seq_record = SeqIO.read(handle, 'fasta')
            fasta_dict[np_num] = str(seq_record.seq)
        return fasta_dict

    def run(self):
        self.check_infiles()
        self.parse_config()
        self.drop_all_tables()
        self.create_all_tables()
        self.register_entries()


# def test():
#     user = 'test_user'
#     password = 'test_password'
#     dbname = 'versus_db'
#     DBE = DataBaseEditor(user, password, dbname)
#     DBE.prepare()
#     # DBE.create_gene_table()
#     # DBE.create_mutation_table()
#     # path = '/Users/moon/DePauw/ITAP/VersUS/result/VersUS.tsv'
#     path = '/Users/moon/DePauw/ITAP/VersUS/src/results/VersUS-test.tsv'
#     vus_dictlist = DBE.parse_vus_data(path)
#     print(vus_dictlist)
#     breakpoint()
#     DBE.register_genes(vus_dictlist)
#     # gene_name = 'POLG'
#     # DBE.get_gene_id(gene_name)
#     # ec_subclass = '/Users/moon/DePauw/ITAP/VersUS/ec_subclass.dat'
#     # ec_class = '/Users/moon/DePauw/ITAP/VersUS/ec_class.txt'
#     ec_file = '/Users/moon/DePauw/ITAP/VersUS/eCNumbersHTML_cleaned.txt'
#     # ec_file = '/Users/moon/DePauw/ITAP/VersUS/ec_exceptions.txt'
#     ec_dict = DBE.parse_EC_numbers(ec_file)
#     DBE.create_ec_table()
#     # ec_dict1 = DBE.parse_EC_subclass(ec_subclass)
#     # ec_dict2 = DBE.parse_EC_class(ec_class)
#     DBE.register_ec_numbers(ec_dict)
#     DBE.close()


def main():
    DBE = DataBaseEditor()
    DBE.run()


if __name__ == '__main__':
    # test()
    main()
