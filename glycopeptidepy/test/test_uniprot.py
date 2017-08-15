import unittest
import requests
from glycopeptidepy.io import uniprot


try:
    is_online = True
    response = requests.get("http://www.uniprot.org/")
    response.raise_for_status()
except:
    is_online = False


def skip_not_online(fn):
    if not is_online:
        return unittest.skip("Not able to reach host")(fn)
    else:
        return fn

@skip_not_online
class UniProtClientTest(unittest.TestCase):
    def test_get(self):
        prot = uniprot.get("P13611")
        self.assertIn("P13611", prot.accessions)
