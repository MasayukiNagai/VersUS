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
import Handlers.util as util


class VersUS:

    def __init__(self):
        self.vus_dict = {}
        self.logger = self.setup_logger('versus_logger', 'versus_logger.log')
        self.time_info = ('VersUS{0:%y%m%d%H%M%S}.log').format(datetime.now())
        self.logger.info('Start VersUS!')


    def argument_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '-c', '--config', type=str, required=True, dest='config',
            help='Required; Specify a config file')
        parser.add_argument(
            '--name', '-n', type=str, required=False, dest='name',
            default=self.time_info,
            help='Required; Specify a suffix for the output')
        parser.add_argument(
            '--clinvar', type=str, required=False, dest='clinvar',
            default=None,
            help='Kickstart: Specify a processed clinvar tsv file')
        parser.add_argument(
            '--blast', type=str, required=False, dest='blast',
            default=None,
            help='Kickstart: Specify a blast xml file')
        parser.add_argument(
            '--vep', type=str, required=False, dest='vep',
            default=None,
            help='Kickstart: Specify a vep tsv file')
        parser.add_argument(
            '--cadd', type=str, required=False, dest='cadd',
            default=None,
            help='Kickstart: Specify a CADD vcf.gz file')
        return parser.parse_args()


    def setup_logger(self, name: str, logfile: str):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s %(funcName)s : %(message)s')
        fh.setFormatter(fh_formatter)
        # creates a file handler that logs messages above INFO level
        sh = logging.StreamHandler()
        sh.setLevel(logging.DEBUG)
        sh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s : %(message)s',
            '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        return logger


    def run(self, config, analysis_id,
            pre_clinvar=None, pre_blast=None, pre_vep=None, pre_cadd=None):
        self.logger.info('Start running the process')
        conf_dict, params_dict = util.parse_config(config)

        # check if config has valid items
        util.check_config_params(params_dict)

        clinvar_file = conf_dict['clinvar']
        genes = conf_dict['genes']
        proteomes = conf_dict['proteomes']
        if pre_blast is None:
            blastp = conf_dict['blastp'] if conf_dict['blastp'] != 'None'\
                 else None
            blastdb = conf_dict['blastdb'] if conf_dict['blastp'] != 'None'\
                 else None
        else:
            blastp, blastdb = None, None
        if pre_vep is None:
            vep = conf_dict['vep'] if conf_dict['vep'] != 'None' else None
        else:
            vep = None
        if pre_cadd is None:
            cadd = True if conf_dict['cadd'] == 'True' else False
        else:
            cadd = None
        paths = [clinvar_file, genes, proteomes, blastp, vep,
                 pre_clinvar, pre_blast, pre_vep, pre_cadd]
        for path in paths:
            if path is not None:
                util.checkpath(path)

        # create correpsonding directories
        interim_dir = os.path.abspath(conf_dict['intermediates'])
        util.make_dir(interim_dir)
        outdir = os.path.abspath(conf_dict['outdir'])
        util.make_dir(outdir)

        seqHandler = SeqHandler()
        seqHandler.setup(conf_dict['email'], conf_dict['apikey'])
        gene_dict = seqHandler.readUniprot_GeneId_EC(genes)

        # parse a ClinvarVariation XML file
        if pre_clinvar is None:
            clinvarHandler = ClinVarHandler(clinvar_file)
            vus_dict = clinvarHandler.readClinVarVariationsXML(
                gene_dict.keys())
        else:
            vus_dict = util.read_tsv_to_dict(pre_clinvar)

        # add EC number and Uniprot id
        vus_dict, uids = seqHandler.add_uniprotId_EC(vus_dict, gene_dict)

        fasta_window = int(params_dict['fasta_window'])
        vus_dict, seq_dict = seqHandler.get_seq(
            vus_dict, proteomes, fasta_window)

        ptmHandler = PTMHandler()
        vus_dict = ptmHandler.addPTM2VUSdict(vus_dict, uids, seq_dict)

        header = list(vus_dict.keys())
        interim_output = os.path.join(interim_dir, f'vus_interim-{analysis_id}.tsv')
        util.write_to_tsv(vus_dict, header, interim_output)

        if blastp:
            blast_input_path = os.path.join(interim_dir, 'blast_input.fasta')
            blast_output_path = os.path.join(interim_dir, 'blast_results.xml')
            blastHandler = BLASTHandler(blastp, blastdb)
            evalue = float(params_dict['evalue'])
            vus_dict = blastHandler.run(
                vus_dict, blast_input_path, blast_output_path, evalue)
        elif pre_blast:
            blastHandler = BLASTHandler()
            vus_dict = blastHandler.run_preprocessed(vus_dict, pre_blast)

        if vep:
            vep_input_path = os.path.join(interim_dir, 'vep_input.tsv')
            vep_output_path = os.path.join(interim_dir, 'vep_results.tsv')
            vepHandler = VEPHandler(vep)
            vus_dict = vepHandler.run(
                vus_dict, vep_input_path, vep_output_path)
        elif pre_vep:
            vepHandler = VEPHandler()
            vus_dict = vepHandler.run_preprocessed(vus_dict, pre_vep)

        if cadd:
            cadd_input_file = os.path.join(interim_dir, 'cadd_input.vcf.gz')
            cadd_output_file = os.path.join(interim_dir, 'cadd_scores.tsv.gz')
            caddHandler = CADDHandler()
            vus_dict = caddHandler.run(
                vus_dict, cadd_input_file, cadd_output_file)
        elif pre_cadd:
            caddHandler = CADDHandler()
            vus_dict = caddHandler.run_preprocessed(vus_dict, pre_cadd)

        header = util.format_header(vus_dict)
        outpath = os.path.join(outdir, f'vus_{analysis_id}.tsv')
        util.write_to_tsv(vus_dict, header, outpath)

        self.logger.info('Finish the process!')


    def main(self):
        start = datetime.now()
        args = self.argument_parser()
        config = args.config
        suffix = args.name
        pre_clinvar = args.clinvar
        pre_blast = args.blast
        pre_vep = args.vep
        pre_cadd = args.cadd

        util.checkpath(config)

        self.run(config, suffix, pre_clinvar, pre_blast, pre_vep, pre_cadd)

        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(f'Running VersUS took {c[0]} minutes {c[1]} seconds')


if __name__ == '__main__':
    versus = VersUS()
    versus.main()
