from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from logging import getLogger

class AlphaFoldHandler:

    def __init__(self):
        self.alphafold_url = 'https://alphafold.ebi.ac.uk/entry/'
        self.logger = getLogger('versus_logger').getChild(__name__)

    def setup(self):
        chrome_options = Options()
        chrome_options.add_argument("--headless")
        self.driver = webdriver.Chrome(ChromeDriverManager().install(),
                                       options=chrome_options)

    def close(self):
        self.driver.close()

    def check_is_valid(self, uniprot_id):
        url = f'{self.alphafold_url}{uniprot_id}'
        self.driver.get(url)
        self.driver.implicitly_wait(1)
        assert 'AlphaFold' in self.driver.title
        try:
            self.driver.find_element_by_id('not-found')
            print(f'{uniprot_id}: not found')
            return False
        except:
            return True

    def check_alphafold(self, vus_dict):
        uniprot_ids = set()
        for var in vus_dict.values():
            uniprot_id = var['uniprot_id']
            uniprot_ids.add(uniprot_id)
        id2url = {}
        print(f'{len(uniprot_ids)} uniprot ids are found in the dict.')
        for i, uniprot_id in enumerate(uniprot_ids):
            is_valid = self.check_is_valid(uniprot_id)
            id2url[uniprot_id] = is_valid
            if i % 100 == 0:
                print(f'{i} ids have been checked.')
        for var in vus_dict.values():
            is_valid = id2url[var['uniprot_id']]
            var['alphafold'] = is_valid
        return vus_dict

    def run(self, vus_dict):
        self.setup()
        self.check_alphafold(vus_dict)
        self.close()

    def test1(self):
        self.setup()
        id1 = 'RANDOM'
        id2 = 'P53396'
        is_valid = self.check_is_valid(id1)
        print(f'uniprot id - {id1}: {is_valid}')
        is_valid = self.check_is_valid(id2)
        print(f'uniprot id - {id2}: {is_valid}')

    def test2(self):
        self.setup()
        vus_dict = {1: {'uniprot_id': 'P53396'},
                    2: {'uniprot_id': 'Q99798'},
                    3: {'uniprot_id': 'P78363'},
                    4: {'uniprot_id': 'RANDOM'},
                    5: {'uniprot_id': 'O00763'},
                    6: {'uniprot_id': 'Q9BZC7'},
                    7: {'uniprot_id': 'P53396'},
                    8: {'uniprot_id': 'Q9UDR5'},
                    9: {'uniprot_id': 'PANDA5'},
                    10: {'uniprot_id': 'P16219'}}
        self.check_alphafold(vus_dict)


if __name__ == '__main__':
    alphafold = AlphaFoldHandler()
    alphafold.test2()
