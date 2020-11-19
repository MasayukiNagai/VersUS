from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import urllib.request
import gzip
from datetime import datetime

class CADDHandler:

    def __init__(self, cadd_input, cadd_output):
        self.cadd_input = cadd_input
        self.cadd_output = cadd_output
        self.cadd_dict = {}
        self.cadd_website_url = "https://cadd.gs.washington.edu/score"


    def setUp(self):
        chrome_options = Options()
        # chrome_options.add_argument("--headless")
        # self.driver = webdriver.Chrome(options = chrome_options)
        self.driver = webdriver.Chrome()


    def upload_CADD_input(self):
        driver = self.driver
        driver.get(self.cadd_website_url)
        assert 'CADD' in driver.title
        upload = driver.find_element_by_xpath("//input[@type='file']")
        upload.send_keys(self.cadd_input)
        version = driver.find_element_by_xpath("//select[@name='version']/option[text()='GRCh38-v1.4']")
        version.click()
        submit = driver.find_element_by_xpath("//input[@type='submit']")
        self.driver.implicitly_wait(5)
        submit.click()


    def check_CADD_upload_succeeded(self):
        if 'success' in self.driver.page_source:
            print('Successfully uploaded a file to CADD')
        elif 'fail' in self.driver.page_source:
            raise NameError('Failed to upload a file to CADD')
        else:
            raise NameError('Cannot tell if a file is successfully uploaded or not')


    def donwload_CADD_results(self):
        url_link = self.check_CADD_output_ready()
        urllib.request.urlretrieve(url_link, self.cadd_output)

    
    def check_CADD_output_ready(self):
        is_ready = False
        start = datetime.datetime.now()
        while(not is_ready):
            url_link = self.driver.find_element_by_xpath(".//a[contains(text(), 'here')]").get_attribute('href')
            self.driver.get(url_link)
            if 'recheck' in self.driver.page_source:
                self.driver.implicitly_wait(60)
            elif 'extension' in self.driver.page_source:
                is_ready = True
                break
            else:
                raise NameError('Cannot tell if a file is successfully uploaded or not')
            lap = datetime.datetime.now()
            time_passed = lap - start
            if time_passed.total_seconds() > 2 * 60 * 60:
                raise NameError('Two hours have passed without CADD scores retrieved. Kill the process.')
        end = datetime.datetime.now()
        time_passed = end - start
        c = divmod(time_passed.days * 86400 + time_passed.seconds, 60)
        print(f'Getting CADD output took {c[0]} minutes {c[1]} seconds')
        return url_link


    def close(self):
        self.driver.close()


    def get_CADD_scores(self):
        self.setUp()
        self.upload_CADD_input()
        self.check_CADD_upload_succeeded()
        self.donwload_CADD_results()
        self.close()
    

    def make_tsv_for_CADD(self, vus_dict):
        header_info = ('CHROM', 'POS', 'ID', 'REF', 'ALT')
        header = '#' + '\t'.join(header_info)
        with open(self.cadd_input, 'w') as f:
            f.write(header + '\n')
            for vus_id in range(len(vus_dict)):
                chrom = vus_dict[vus_id]['chr']
                pos = str(vus_dict[vus_id]['start'])
                ref = vus_dict[vus_id]['referenceAllele']
                alt = vus_dict[vus_id]['alternateAllele']
                info_tup = (chrom, pos, str(vus_id), ref, alt)
                info = '\t'.join(info_tup)
                f.write(info + '\n')
    

    def read_CADD_results(self):
        with gzip.open(self.cadd_output, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                ls = line.split('\t')
                chrom = ls[0].strip()
                pos = ls[1].strip()
                ref = ls[2].strip()
                alt = ls[3].strip()
                phred = ls[5].strip()
                key = ref + pos + alt
                if chrom not in self.cadd_dict.keys():
                    self.cadd_dict[chrom] = {}
                self.cadd_dict[chrom][key] = phred
        length = 0
        for chrom in self.cadd_dict.keys():
            length += len(self.cadd_dict[chrom].keys())
        print(f'Length of CADD dict: {length}')
        return self.cadd_dict


    def add_cadd_results(self, vus_dict: dict):
        unfound_cadd = {}
        for vus_id in vus_dict:
            chrom = vus_dict[vus_id]['chr']
            ref = vus_dict[vus_id]['referenceAllele']
            pos = vus_dict[vus_id]['start']
            alt = vus_dict[vus_id]['alternateAllele']
            key = ref + pos + alt
            try:
                cadd_score = self.cadd_dict[chrom][key]
            except:
                c_acc = vus_dict[vus_id]['ClinVar_accession']
                unfound_cadd[c_acc] = key
                cadd_score = None
            vus_dict[vus_id]['CADD_score'] = cadd_score
        print('Unfound cadd: ', unfound_cadd)
        return vus_dict

    
    def run(self, vus_dict):
        self.make_tsv_for_CADD(vus_dict)
        self.get_CADD_scores()
        self.read_CADD_results()
        vus_dict = self.add_cadd_results(vus_dict)
        return vus_dict
        

# cadd_input = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/CADD_sample_input.vcf'
# cadd_output = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/CADD_sample_scores.tsv.gz'
# ch = CADDHandler(cadd_input, cadd_output)
# ch.get_CADD_scores()

# need to make a method to change a relative path to abs path