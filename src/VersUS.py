import os
import logging
import argparse
from datetime import datetime
from Handlers.SeqHandler import *
from Handlers.ClinVarHandler import *
from Handlers.BLASTHandler import *
from Handlers.CADDHandler import *
from Handlers.VEPHandler import *
from Handlers.FileHandler import *


class VersUS:

    def __init__(self):
        self.vus_dict = {}
        self.logger = self.setup_logger('versus_logger', 'versus_logger.log')
        time_info = ('VersUS{0:%y%m%d%H%M%S}.log').format(datetime.now())
        self.logger.info('Start VersUS!')


    def argument_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config', '-c', nargs=1, type=str, required=True,
                            help='Required; Specify a config file.')
        args = parser.parse_args()
        return args


    def setup_logger(self, name: str, logfile: str):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level 
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter('[%(asctime)s] %(levelname)s %(filename)s %(funcName)s : %(message)s')
        fh.setFormatter(fh_formatter)
        # creates a file handler that logs messages above INFO level 
        sh = logging.StreamHandler()
        sh.setLevel(logging.DEBUG)
        sh_formatter = logging.Formatter('[%(asctime)s] %(levelname)s %(filename)s : %(message)s', '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        return logger


    def run(self, config):
        self.logger.info('Start running the process')
        general_dict, params_dict = parse_config(config)

        # check if config has valid items
        check_config_general(general_dict)
        check_config_params(params_dict)

        genes = general_dict['genes']
        proteomes = general_dict['proteomes']
        clinvar_variations = general_dict['variations']
        blast = os.path.abspath(general_dict['blast']) if general_dict['blast'] != 'None' else None
        vep = os.path.abspath(general_dict['vep']) if general_dict['vep'] != 'None' else None
        cadd = True if general_dict['cadd'] == 'True' else False

        # create correpsonding directories
        intermediates_dir = os.path.abspath(general_dict['intermediates'])
        make_dir(intermediates_dir)
        outdir = os.path.abspath(general_dict['outdir'])
        make_dir(outdir)

        seqHandler = SeqHandler(genes, proteomes)
        genes_dict = seqHandler.readHumanGenesEC()

        clinvarHandler = ClinVarHandler(clinvar_variations)
        vus_dict = clinvarHandler.readClinVarVariationsXML(genes_dict)
        
        fasta_window = int(params_dict['fasta_window'])
        vus_dict = seqHandler.get_seq(vus_dict, fasta_window)
        
        if blast:
            blast_input_path = os.path.join(intermediates_dir, 'blast_vus_input.fasta')
            blast_output_path = os.path.join(intermediates_dir, 'blast_vus_results.xml')
            blastHandler = BLASTHandler(blast, blast_input_path, blast_output_path)
            evalue = float(params_dict['evalue'])
            vus_dict = blastHandler.run(vus_dict, evalue)
            
            header = format_header(vus_dict)
            outpath = os.path.join(intermediates_dir, 'vus_blast.tsv')
            write_to_tsv(vus_dict, header, outpath)
        
        if vep:
            vep_input_path = os.path.join(intermediates_dir, 'vep_vus_input.tsv')
            vep_output_path = os.path.join(intermediates_dir, 'vep_vus_results.tsv')
            vepHandler = VEPHandler(vep, vep_input_path, vep_output_path)
            vus_dict = vepHandler.run(vus_dict)

            header = format_header(vus_dict)
            outpath = os.path.join(intermediates_dir, 'vus_vep.tsv')
            write_to_tsv(vus_dict, header, outpath)
        
        if cadd:
            cadd_input_file = os.path.join(intermediates_dir, 'cadd_vus_input.tsv')
            cadd_output_file = os.path.join(intermediates_dir, 'cadd_vus_scores.tsv.gz')
            caddHandler = CADDHandler(cadd_input_file, cadd_output_file)
            vus_dict = caddHandler.run(vus_dict)

        header = format_header(vus_dict)
        outpath = os.path.join(outdir, 'VersUS.tsv')
        write_to_tsv(vus_dict, header, outpath)

        self.logger.info('Finish the process!')


    def main(self):
        start = datetime.now()
        args = self.argument_parser()
        config = args.config[0]

        checkpath(config)

        self.run(config)

        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Running VersUS took {c[0]} minutes {c[1]} seconds')
        
        
if __name__ == '__main__':
    versus = VersUS()
    versus.main()