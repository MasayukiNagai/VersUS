import os
import sys
import logging


class VersUS:

    def __init__(self):
        pass


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


    def make_output_dir(self, outdir: str):
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except OSError:
                sys.exit('Cannot make an output directory.')
        else:
            pass


    def run(self, logger, outdir: str):
        pass


    def main(self):
        pass
