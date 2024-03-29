import os
from datetime import datetime
from logging import getLogger


class VEPHandler:

    def __init__(self, vep='vep'):
        self.vep = vep
        self.logger = getLogger('versus_logger').getChild(__name__)

    def make_ordered_vus_dict_for_vep(self, vus_dict):
        ordered_dict = {}
        ct_none = 0
        self.logger.debug(f'Lenght of original_dict: {len(vus_dict)}')
        for vus_id in vus_dict.keys():
            chrom = vus_dict[vus_id]['chr']
            try:
                start = int(vus_dict[vus_id]['start'])
                stop = int(vus_dict[vus_id]['stop'])
            except Exception:
                ct_none += 1
                continue
            ref = vus_dict[vus_id]['referenceAllele']
            alt = vus_dict[vus_id]['alternateAllele']
            if ref is None or alt is None:
                ct_none += 1
                continue
            if chrom not in ordered_dict.keys():
                ordered_dict[chrom] = {}
            if start not in ordered_dict[chrom].keys():
                ordered_dict[chrom][start] = [{'ref': ref, 'alt': alt}]
            else:
                ordered_dict[chrom][start].append({'ref': ref, 'alt': alt})
            if start != stop:
                print(f'Start != Stop: {start}-{stop}')
        ct = 0
        for chrom in ordered_dict.keys():
            for pos in ordered_dict[chrom].keys():
                ct += len(ordered_dict[chrom][pos])
        self.logger.debug(f'Length of ordered_dict: {ct}')
        self.logger.debug(f'Number of vus with any info missing: {ct_none}')
        return ordered_dict

    def make_tsv_ordered_for_vep(self, vus_ordered_dict, vep_infile):
        chrom_order = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                       '11', '12', '13', '14', '15', '16', '17', '18', '19',
                       '20', '21', '22', 'X')
        with open(vep_infile, 'w') as f:
            for chrom in chrom_order:
                if chrom not in vus_ordered_dict.keys():
                    continue
                for pos in sorted(vus_ordered_dict[chrom].keys()):
                    for mut in vus_ordered_dict[chrom][pos]:
                        info = (
                            chrom, str(pos), str(pos), mut['ref'], mut['alt'])
                        f.write('\t'.join(info) + '\n')

    def make_tsv_for_vep(self, vus_dict, vep_infile):
        with open(vep_infile, 'w') as f:
            for vus_id in range(0, len(vus_dict)):
                info = (
                    vus_dict[vus_id]['chr'],
                    vus_dict[vus_id]['start'],
                    vus_dict[vus_id]['stop'],
                    vus_dict[vus_id]['referenceAllele'],
                    vus_dict[vus_id]['alternateAllele'])
                f.write('\t'.join(info) + '\n')

    def vep_locally(self, infile, outfile):
        self.logger.info('Start running VEP')
        start = datetime.now()
        cmd = f'{self.vep} '\
            + f'-i {infile} '\
            + f'-o {outfile} '\
            + '--cache ' + '--af_gnomadg '\
            + '--no_check_variants_order '\
            + '--flag_pick ' + '--tab '\
            + '--force_overwrite ' + '--offline'
            # + '--appris ' + '--canonical '\
            # + '--sift p ' + '--polyphen p '\
        vep_cmd = os.system(cmd)
        end = datetime.now()
        time = end - start
        c = divmod(time.days * 86400 + time.seconds, 60)
        self.logger.info(cmd + ' : ran with exit code %d' %vep_cmd)
        self.logger.info(f'Running Vep took {c[0]} minutes {c[1]} seconds')

    # def crop_vep_output(self, vep_file_path: str, outfile_path: str):
    #     cmd = "cat " + vep_file_path + " | "\
    #         + "grep -v '##'" + " | "\
    #         + "sed 's/^#\(.*\)/\\1/'" +  " >| "\
    #         + outfile_path
    #     crop_cmd = os.system(cmd)
    #     print(cmd + ' : ran with exit code %d' %crop_cmd)

    def read_vep_output(self, outfile):
        vep_dict = {}
        with open(outfile, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('##'):
                    continue
                elif line.startswith('#'):
                    header = line.split('\t')
                    location_i = header.index('Location')
                    alt_i = header.index('Allele')
                    pick_i = header.index('PICK')
                    gnomadAF_i = header.index('gnomADg_AF')
                    self.logger.debug(f'VEP output header: {header}')
                    continue
                data = line.split('\t')
                if len(header) != len(data):
                    self.logger.warning(
                        'Error: header and data have different line length')
                    continue
                if data[pick_i] != '1':
                    continue
                chrom, pos = data[location_i].split(':')
                alt = data[alt_i]
                try:
                    gnomadAF = float(data[gnomadAF_i])
                except Exception:
                    gnomadAF = None
                if chrom not in vep_dict.keys():
                    vep_dict[chrom] = {}
                key = pos + '>' + alt
                vep_dict[chrom][key] = gnomadAF
        ct = 0
        for chrom in vep_dict:
            ct += len(vep_dict[chrom])
        self.logger.debug(f'Number of items in the vep dict: {ct}')
        return vep_dict

    def add_vep_output(self, vus_dict, vep_dict):
        unfound_af = {}
        for vus in vus_dict.values():
            chrom = vus['chr']
            try:
                key = str(vus['start']) + '>' + vus['alternateAllele']
                gnomAD_AF = vep_dict[chrom][key]
            except Exception:
                gnomAD_AF = None
                c_acc = vus['ClinVar_accession']
                unfound_af[c_acc] = chrom + ':' + key
            vus['gnomAD_AF'] = gnomAD_AF
        self.logger.debug(f'Unfound gnomAD_AF: {unfound_af}')
        return vus_dict

    def run(self, vus_dict, vep_infile, vep_outfile):
        vus_ordered_dict = self.make_ordered_vus_dict_for_vep(vus_dict)
        self.make_tsv_ordered_for_vep(vus_ordered_dict, vep_infile)
        self.vep_locally(vep_infile, vep_outfile)
        vep_dict = self.read_vep_output(vep_outfile)
        vus_dict = self.add_vep_output(vus_dict, vep_dict)
        return vus_dict

    def run_preprocessed(self, vus_dict, vep_outfile):
        vep_dict = self.read_vep_output(vep_outfile)
        vus_dict = self.add_vep_output(vus_dict, vep_dict)
        return vus_dict
