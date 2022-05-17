import argparse
import logging
from Handlers.SeqHandler import *
from Handlers.ClinVarHandler import *


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genelist', '-g', type=str, required=True,
                        dest='gene', help='Required; uniprot gene file.')
    parser.add_argument('--xml', '-x', type=str, required=True,
                        dest='xml', help='Required; ClinVar xml file.')
    args = parser.parse_args()
    return args


def setup_logger(name: str, logfile: str):
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


def main():
    args = argument_parser()
    genelist = args.gene
    xml = args.xml
    logger = setup_logger('versus_logger', 'vus_stats.log')
    seqHandler = SeqHandler(genelist, None)
    genes_dict = seqHandler.readUniprot_GeneId_EC()
    clinvarHandler = ClinVarHandler(xml)
    vus_dict = clinvarHandler.readClinVarVariationsXML(genes_dict.keys())


if __name__ == '__main__':
    main()
