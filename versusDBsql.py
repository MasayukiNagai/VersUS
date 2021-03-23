import mysql.connector

class DataBaseEditor:
    def __init__(self, user, passwd, dbname):
        self.host = 'localhost'
        self.user = user
        self.passwd = passwd
        self.dbname = dbname
        self.gene_table = 'Gene'
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
        name_type_dict = {'gene_id': 'INT PRIMARY KEY',
                          'gene_name_short': 'VARCHAR(10) ' + nn,
                          'gene_name_full': 'VARCHAR(255) ' + nn,
                          'chrom': 'VARCHAR(20) ' + nn,
                          'EC_number': 'VARCHAR(20) ' + nn,
                          'EC_id': 'VARCHAR(20) ' + nn}
        self.create_table(self.gene_table, name_type_dict)
    

    def create_mutation_table(self):
        nn = 'NOT NULL'
        name_type_dict = {'mutation_id': 'INT PRIMARY KEY',
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


