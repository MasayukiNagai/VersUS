import argparse
from datetime import datetime
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


def main():
    args = argument_parser()
    genelist = args.gene
    xml = args.xml
    seqHandler = SeqHandler(genelist, None)
    genes_dict = seqHandler.readUniprot_GeneId_EC()
    clinvarHandler = ClinVarHandler(xml)
    vus_dict = clinvarHandler.readClinVarVariationsXML(genes_dict.keys())


if __name__ == '__main__':
    main()
