import os
import logging
import argparse
from datetime import datetime
from Handlers.SeqHandler import SeqHandler
from Handlers.ClinVarHandler import ClinVarHandler
from Handlers.BLASTHandler import BLASTHandler
from Handlers.CADDHandler import CADDHandler
from Handlers.VEPHandler import VEPHandler
from Handlers.PTMHandler import PTMHandler
from Handlers.FileHandler import *


class VersUS:

    def __init__(self):
        self.vus_dict = {}
        self.logger = self.setup_logger('versus_logger', 'versus_logger.log')
        time_info = ('VersUS{0:%y%m%d%H%M%S}.log').format(datetime.now())
        self.logger.info('Start VersUS!')


    def argument_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config', '-c', type=str, required=True,
                            dest='config', help='Required; Specify a config file.')
        parser.add_argument('--name', '-n', type=str, required=True,
                            dest='name', help='Required; Specify analysis-ID that is added to the end of outputs.')
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


    def run(self, config, analysis_id):
        self.logger.info('Start running the process')
        conf_dict, params_dict = parse_config(config)

        # check if config has valid items
        check_config_params(params_dict)

        clinvar_file = conf_dict['clinvar']
        genes = conf_dict['genes']
        proteomes = conf_dict['proteomes']
        blast = os.path.abspath(conf_dict['blast']) if conf_dict['blast'] != 'None' else None
        vep = os.path.abspath(conf_dict['vep']) if conf_dict['vep'] != 'None' else None
        cadd = True if conf_dict['cadd'] == 'True' else False
        for path in [clinvar_file, genes, proteomes, blast, vep]:
            if path is not None:
                checkpath(path)

        # create correpsonding directories
        intermediates_dir = os.path.abspath(conf_dict['intermediates'])
        make_dir(intermediates_dir)
        outdir = os.path.abspath(conf_dict['outdir'])
        make_dir(outdir)

        seqHandler = SeqHandler(genes, proteomes)
        gene_dict = seqHandler.readUniprot_GeneId_EC()

        # parse a ClinvarVariation XML file
        clinvarHandler = ClinVarHandler(clinvar_file)
        vus_dict = clinvarHandler.readClinVarVariationsXML(gene_dict.keys())

        # add EC number and Uniprot id
        vus_dict = seqHandler.add_uniprotId_EC(vus_dict)

        fasta_window = int(params_dict['fasta_window'])
        vus_dict = seqHandler.get_seq(vus_dict, fasta_window)

        ptmHandler = PTMHandler()
        vus_dict = ptmHandler.addPTM2VUSdict(vus_dict, gene_dict)

        # header = format_header(vus_dict)
        # intermediate_output = os.path.join(intermediates_dir, f'vus_intermediate-{analysis_id}.tsv')
        # write_to_tsv(vus_dict, header, intermediate_output)

        if blast:
            blast_input_path = os.path.join(intermediates_dir, 'blast_input.fasta')
            blast_output_path = os.path.join(intermediates_dir, 'blast_results.xml')
            blastHandler = BLASTHandler(blast, blast_input_path, blast_output_path)
            evalue = float(params_dict['evalue'])
            vus_dict = blastHandler.run(vus_dict, evalue)
            # header = format_header(vus_dict)
            # write_to_tsv(vus_dict, header, intermediate_output)

        if vep:
            vep_input_path = os.path.join(intermediates_dir, 'vep_input.tsv')
            vep_output_path = os.path.join(intermediates_dir, 'vep_results.tsv')
            vepHandler = VEPHandler(vep, vep_input_path, vep_output_path)
            vus_dict = vepHandler.run(vus_dict)
            # header = format_header(vus_dict)
            # write_to_tsv(vus_dict, header, intermediate_output)

        if cadd:
            cadd_input_file = os.path.join(intermediates_dir, 'cadd_input.vcf.gz')
            cadd_output_file = os.path.join(intermediates_dir, 'cadd_scores.tsv.gz')
            caddHandler = CADDHandler(cadd_input_file, cadd_output_file)
            vus_dict = caddHandler.run(vus_dict)

        header = format_header(vus_dict)
        outpath = os.path.join(outdir, f'vus_{analysis_id}.tsv')
        write_to_tsv(vus_dict, header, outpath)

        self.logger.info('Finish the process!')


    def main(self):
        start = datetime.now()
        args = self.argument_parser()
        config = args.config
        analysis_id = args.name

        checkpath(config)

        self.run(config, analysis_id)

        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Running VersUS took {c[0]} minutes {c[1]} seconds')



if __name__ == '__main__':
    versus = VersUS()
    versus.main()
