import argparse
import gzip
from lxml import etree


# need to change global WFILE stuff later
class VariationHandlerSpecific(object):

    def __init__(self, accession, outfile=None):
        self.is_accession = False
        self.accession = accession
        self.ct = 0
        print(f"Accession: {accession}")
        self.wfile = outfile if outfile else f'ClinVar_{accession}.xml'
        open(self.wfile, "w").close()  # erace contents

    def start(self, tag, attrs):
        # global WFILE
        if (tag == 'VariationArchive') \
            and (attrs.get('Accession') == self.accession):
            self.is_accession = True
            print('The specific variation is found: ' + str(self.ct))
        if self.is_accession:
            if len(attrs.keys()) == 0:
                with open(self.wfile, 'a') as f:
                    f.write('<' + tag)
            else:
                for i, t in enumerate(attrs.keys()):
                    if i == 0:
                        with open(self.wfile, 'a') as f:
                            f.write('<' + tag + ' ')
                    elif i != len(attrs.keys()) - 1:
                        with open(self.wfile, 'a') as f:
                            f.write(t + '="' + attrs.get(t) + '"' + ' ')
                    else:
                        with open(self.wfile, 'a') as f:
                            f.write(t + '="' + attrs.get(t) + '"')
            with open(self.wfile, 'a') as f:
                f.write('>')

    def end(self, tag):
        # global WFILE
        if self.is_accession:
            with open(self.wfile, 'a') as f:
                f.write('</' + tag + '>')
        if tag == 'VariationArchive':
            self.ct += 1
            if self.ct % 100000 == 0:
                print(self.ct)
        if self.is_accession and tag == 'VariationArchive':
            self.is_accession = False
            # WFILE.close()
            print('The subnode file is completed')


    def data(self, data):
        # global WFILE
        if data is not None:
            if self.is_accession and data != "":
                with open(self.wfile, 'a') as f:
                    f.write(data)

    def close(self):
        # WFILE.close()
        print('The xml file is closed')


# read xml file of variations from ClinVar
# return dataframe and write to a csv file
def readClinVarVariationsXMLSpecific(clinvar_xml, accession, outfile):
    parser = etree.XMLParser(target=VariationHandlerSpecific(accession, outfile))
    with gzip.open(clinvar_xml, 'rt') as handle:
        etree.parse(handle, parser)


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--accession', type=str, required=True, dest='acc',
        help='Required; specify an accession number')
    parser.add_argument(
        '-x', '--xml', type=str, required=True, dest='xml',
        help='Required; specify an xml file')
    parser.add_argument(
        '-o', '--out', type=str, required=False, dest='out',
        default=None,
        help=('Optional; specify an output path. '
              'Default "ClinVar_{accession}.xml"'))
    return parser.parse_args()


if __name__ == '__main__':
    args = argument_parser()
    xml = args.xml
    acc = args.acc
    outfile = args.out
    readClinVarVariationsXMLSpecific(xml, acc, outfile)
