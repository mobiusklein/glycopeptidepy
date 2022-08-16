import io

try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen

from . import fasta

peff_api_url = 'https://api.nextprot.org/export/entry/NX_%s.peff'


def download_peff(accession):
    url = peff_api_url % accession
    return fasta.PEFFReader(io.BytesIO(urlopen(url).read()), index=True)
