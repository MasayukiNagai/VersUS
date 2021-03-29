import mysql.connector

aaMapOneToThree = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 
                   'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 
                   'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 
                   'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}

cinvar_link_variation = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/'

class DataBaseEditor:
    def __init__(self, user, passwd, dbname):
        self.host = 'localhost'
        self.user = user
        self.passwd = passwd
        self.dbname = dbname
        self.gene_table = 'Gene2'
        self.mutation_table = 'Mutation'
        # self.clinical_significance = 'Clinical_Significane'
        # self.pdb = 'PDB'
        self.gene_id = 0
        self.mutation_id = 0

    
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
        nn = 'NOT NULL'
        name_type_dict = {'gene_id': 'SERIAL PRIMARY KEY',
                          'gene_name_short': 'VARCHAR(10) ' + nn,
                          'gene_name_full': 'VARCHAR(255) ' + nn,
                          'chrom': 'VARCHAR(20) ' + nn,
                          'EC_number': 'VARCHAR(20) ' + nn}
                        #   'EC_id': 'VARCHAR(20) ' + nn
        self.create_table(self.gene_table, name_type_dict)
    

    def create_mutation_table(self):
        nn = 'NOT NULL'
        name_type_dict = {'mutation_id': 'SERIAL PRIMARY KEY',
                          'gene_id': 'INT ' + nn,
                          'ref_pos_alt': 'VARCHAR(20)' + nn,
                          'clinvar_link': 'VARCHAR(255)' + nn,
                          'clinical_significance': 'VARCHAR(255)' + nn,
                          'CADD_score': 'FLOAT',
                          'gnomAD_AF': 'FLOAT',
                          'pdb': 'VARCHAR(20)'}
        self.create_table(self.mutation_table, name_type_dict)
    

    def create_all_tables(self):
        # connect to mysql database
        self.prepare()
        # create tables
        self.create_gene_table()
        self.create_mutation_table()
        self.close()
        # close connection
        self.close()

    
    def insert_items(self, table, keytup, values_tuplist):
        query1 = 'INSERT INTO ' + table + ' '
        body1 = '(' + ', '.join(keytup) + ') VALUES '
        values_query_list = []
        for values in values_tuplist:
            val_q = '(' + ', '.join(values) + ')'
            values_query_list.append(val_q)
        body2 = ', '.join(values_query_list) + ';'
        query = query1 + body1 + body2
        self.cur.execute(query)
        self.cnx.commit()
        print(query)
    

    def str_editor(self, strings):
        if strings == 'NULL' or strings == 'null':
            return strings
        if type(strings) is str:
            return "'" + strings + "'"
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


    def register_vus(self, vus_tsv):
        vus_dictlist = self.parse_vus_data(vus_tsv)
        self.register_genes(vus_dictlist)
        self.register_mutations(vus_dictlist)


    def register_genes(self, vus_dictlist):
        gene_dict = dict()
        for vus_dict in vus_dictlist:
            if vus_dict['gene_id'] not in gene_dict.keys():
                gene_dict[vus_dict['gene_id']] = {'gene_name_short': vus_dict['gene_id'],
                                                  'gene_name_full': vus_dict['gene_name'],
                                                  'chrom': vus_dict['chr'],
                                                  'EC_number': vus_dict['EC_number']}
        keytup = ('gene_name_short', 'gene_name_full', 'chrom', 'EC_number')
        gene_values = self.get_tuplist_for_gene(keytup, gene_dict)
        self.insert_items(self.gene_table, keytup, gene_values)


    def register_mutations(self, vus_dictlist):
        mutation_dict = dict()
        for vus_dict in vus_dictlist:
            if vus_dict['gene_id'] not in mutation_dict.keys():
                mutation_dict[vus_dict['gene_id']] = []
            mutation_dict[vus_dict['gene_id']].append(self.make_mut_dict(vus_dict))
        keytup = ('gene_id', 'ref_pos_alt', 'clinvar_link', 'clinical_significance', 'CADD_score', 'gnomAD_AF', 'pdb')
        vus_values = self.get_tuplist_for_mut(keytup, mutation_dict)
        self.insert_items(self.mutation_table, keytup, vus_values)

    
    def make_mut_dict(self, vus_dict):
        gene_id = self.get_gene_id(vus_dict['gene_id'])
        ref_pos_alt = vus_dict['missense_variation']
        # ref_pos_alt = aaMapOneToThree(vus_dict['ref']) + vus_dict['pos'] +  aaMapOneToThree(vus_dict['alt'])
        clinvar_link = cinvar_link_variation + vus_dict['ClinVar_accession']
        mut = {'gene_id': gene_id,
               'ref_pos_alt': ref_pos_alt, 
               'clinvar_link': clinvar_link,
               'clinical_significance': vus_dict['clinical_significance'],
               'CADD_score': vus_dict['CADD_score'],
               'gnomAD_AF': vus_dict['gnomAD_AF'],
               'pdb': vus_dict['pdb_ID']}
        return mut


    def get_gene_id(self, gene_name):
        query = 'SELECT gene_id FROM ' + self.gene_table\
               + ' WHERE gene_name_short = ' + self.str_editor(gene_name)
        self.cur.execute(query)
        gene_id_list = self.cur.fetchall()
        assert len(gene_id_list) == 1, f'"{gene_name}" does not exist in Gene table'
        gene_id = gene_id_list[0][0]
        return gene_id


    def parse_vus_data(self, vus_tsv):
        with open(vus_tsv, 'r') as f:
            header = f.readline().rstrip()
            keys = header.split('\t')
            keytup = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'FASTA_window', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
            assert len(keys) == len(keytup)
            vus_dictlist = []
            for line in f:
                data = line.rstrip().split('\t')
                assert len(data) == len(keytup)
                vus_dict = dict(zip(keytup, data))
                for key in keytup:
                    if vus_dict[key] == 'None':
                        vus_dict[key] = 'NULL'
                vus_dictlist.append(vus_dict)
        return vus_dictlist


def main():
    user = 'test_user'
    password = 'test_password'
    dbname = 'versus_db'
    DBE = DataBaseEditor(user, password, dbname)
    DBE.prepare()
    DBE.create_gene_table()
    DBE.create_mutation_table()
    path = '/Users/moon/DePauw/ITAP/VersUS/result/VersUS.tsv'
    DBE.register_vus(path)
    # vus_dictlist = DBE.parse_vus_data(path)
    # DBE.register_genes(vus_dictlist)
    # gene_name = 'POLG'
    # DBE.get_gene_id(gene_name)
    DBE.close()


if __name__ == '__main__':
    main()
