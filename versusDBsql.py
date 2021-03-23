import mysql.connector

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
                          'ref': 'CHAR(1) ' + nn,
                          'alt': 'CHAR(1) ' + nn,
                          'pos': 'INT ' + nn,
                          'clinical_significance_id': 'INT ' + nn,
                          'NP_accession': 'VARCHAR(255)',
                          'Clinvar_accession': 'VARCHAR(255) ' + nn,
                          'vcf_ref': 'CHAR(1) ' + nn,
                          'vcf_alt': 'CHAR(1) ' + nn,
                          'vcf_pos': 'INT ' + nn,
                          'gnomAD_AF': 'FLOAT',
                          'CADD_score': 'FLOAT',
                          'BLAST_evalue': 'FLOAT',
                          'pdb_id': 'INT',
                          'hit_from': 'INT',
                          'hit_to': 'INT',
                          'FASTA_window': 'VARCHAR(255)'}
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
        print(query)
    

    def str_editor(self, strings):
        if strings == 'NULL' or strings == 'null':
            return strings
        if type(strings) is str:
            return "'" + strings + "'"
        else:
            return str(strings)
    

    # def get_tuplist_for_vus(self, keytup, vus_dictlist):
    #     tuplist = []
    #     for vus_dict in vus_dictlist:
    #         tup = tuple([vus_dict[key] for key in keytup])
    #         tuplist.append(tup)
    #     return tuplist


    def get_tuplist_for_gene(self, keytup, gene_dict):
        tuplist = []
        for gene_id in gene_dict.keys():
            tup = tuple([self.str_editor(gene_dict[gene_id][key]) for key in keytup])
            tuplist.append(tup)   
        return tuplist


    def register_vus(self, vus_tsv):
        vus_dictlist = self.parse_vus_data(vus_tsv)


    # def register_mutations(self):
    #     # Keys for items to be registered (excluded FASTA window)
    #     keytup = ('gene_id', 'gene_name', 'clinical_significance', 'EC_number', 'missense_variation', 'NP_accession', 'ClinVar_accession', 'gnomAD_AF', 'CADD_score', 'chr', 'start', 'stop', 'referenceAllele', 'alternateAllele', 'pdb_ID', 'BLAST_evalue', 'hit_from', 'hit_to')
    #     vus_values = self.get_tuplist_for_vus(keytup, vus_dictlist)
    #     self.insert_items(self.mutation_table)

    
    def register_genes(self, vus_dictlist):
        gene_dict = dict()
        for vus_dict in vus_dictlist:
            if vus_dict['gene_id'] not in gene_dict.keys():
                gene_dict[vus_dict['gene_id']] = {}
                gene_dict[vus_dict['gene_id']]['gene_name_short'] = vus_dict['gene_id']
                gene_dict[vus_dict['gene_id']]['gene_name_full'] = vus_dict['gene_name']
                gene_dict[vus_dict['gene_id']]['chrom'] = vus_dict['chr']
                gene_dict[vus_dict['gene_id']]['EC_number'] = vus_dict['EC_number']
        keytup = ('gene_name_short', 'gene_name_full', 'chrom', 'EC_number')
        gene_values = self.get_tuplist_for_gene(keytup, gene_dict)
        self.insert_items(self.gene_table, keytup, gene_values)
  
    
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
    path = '/Users/moon/DePauw/ITAP/VersUS/result/VersUS.tsv'
    vus_dictlist = DBE.parse_vus_data(path)
    # DBE.register_genes(vus_dictlist)
    DBE.close()


if __name__ == '__main__':
    main()
