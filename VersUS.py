import os
import sys
import logging


class VersUS:

    def __init__(self):
        self.vus_dict = {}
        pass


    def argument_parser(self):
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
    

    def write_to_tsv(self, vus_dict: dict, header: tuple, outfile_path: str):
        with open(outfile_path, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for vus_id in vus_dict:
                info = []
                for item in header:
                    info.append(str(vus_dict[vus_id][item]))
                f.write('\t'.join(info) + '\n')  


    def write_to_csv(self, vus_dict: dict, header: tuple, outfile_path: str):
        with open(outfile_path, 'w') as f:
            f.write(','.join(header) + '\n')
            for vus_id in vus_dict:
                info = []
                for item in header:
                    info.append(str(vus_dict[vus_id][item]))
                f.write(','.join(info) + '\n')

    
    def read_tsv_to_dict(self, tsv_path: str):
        vus_dict = {}
        with open(tsv_path, 'r') as f:
            header = f.readline().split('\t')
            header = [item.strip() for item in header]
            print(header)
            for i, line in enumerate(f):
                vus_dict[i] = {}
                ls = line.split('\t')
                if (len(header) != len(ls)):
                    print(f'The length is different: {line}')
                    continue
                for j, item in enumerate(header):
                    vus_dict[i][item] = ls[j].strip()
        return vus_dict 


    def read_csv_to_dict(self, csv_path: str):
        vus_dict = {}
        with open(csv_path, 'r') as f:
            header = f.readline().split(',')
            header = [item.strip() for item in header]
            for i, line in enumerate(f):
                vus_dict[i] = {}
                ls = line.split(',')
                for j, item in enumerate(header):
                    vus_dict[i][item] = ls[j].strip()
        return vus_dict


    def run(self, logger, outdir: str):
        pass


    def main(self):
        pass
