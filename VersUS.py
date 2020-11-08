import os
import logging
import argparse
from datetime import datetime
from ClinVarHandler import *
from BLASTHandler import *
from CADDHandler import *
from VEPHandler import *
from FileHandler import *


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
        sh.setLevel(logging.INFO)
        sh_formatter = logging.Formatter('[%(asctime)s] %(levelname)s : %(message)s', '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        return logger


    def run(self, config):
        self.logger.info('Start the process')


    def main(self):
        args = self.argument_parser()
        config = args.config[0]
        
