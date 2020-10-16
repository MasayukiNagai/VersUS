from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import urllib.request
import datetime

class WebsiteHandler:

    def setUp(self):
        chrome_options = Options()
        chrome_options.add_argument("--headless")
        self.driver = webdriver.Chrome(options = chrome_options)
        self.driver = webdriver.Chrome()
    
    def get_CADD_scores(self, cadd_inputfile: str, cadd_outfile: str):
        self.setUp()
        self.upload_CADD_input(cadd_inputfile)
        self.check_CADD_upload_succeeded()
        self.donwload_CADD_results(cadd_outfile)
        self.close()

    def upload_CADD_input(self, cadd_inputfile: str):
        driver = self.driver
        driver.get("https://cadd.gs.washington.edu/score")
        assert 'CADD' in driver.title
        upload = driver.find_element_by_xpath("//input[@type='file']")
        upload.send_keys(cadd_inputfile)
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

    def donwload_CADD_results(self, outfile_path):
        url_link = self.check_CADD_output_ready()
        urllib.request.urlretrieve(url_link, outfile_path)
    
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
            if time_passed.total_seconds() > 60 * 60:
                raise NameError('One hour has passed')
        end = datetime.datetime.now()
        time_passed = end - start
        c = divmod(time_passed.days * 86400 + time_passed.seconds, 60)
        print(f'Getting CADD output took {c[0]} minutes {c[1]} seconds')
        return url_link

    def close(self):
        self.driver.close()


# cadd_inputfile = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/CADD_sample_input.vcf'
# cadd_outfile = '/Users/moon/DePauw/ITAP/ClinvarSorting/data/CADD/test_results.tsv.gz'
# wh = WebsiteHandler()
# wh.get_CADD_scores(cadd_inputfile, cadd_outfile)

# need to make a method to change a relative path to abs path