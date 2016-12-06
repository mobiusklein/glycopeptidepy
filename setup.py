import re
from setuptools import setup, find_packages, Extension


def version():
    return re.sub(r'[\s\'"\n]', '', open("glycopeptidepy/version.py").readline().split("=")[1])


required = []
with open('requirements.txt') as f:
    required = f.read().splitlines()


def run_setup(include_cext=True):
    setup(
        name='glycopeptidepy',
        version=version(),
        packages=find_packages(),
        install_requires=required,
        include_package_data=True,
        zip_safe=False,
        package_data={
            'glycopeptidepy': ["*.csv", "*.xml", "*.json", "data/*.csv"],
            'glycopeptidepy.structure': ["structure/data/*.csv", "structure/data/*.json"]
        },
    )


run_setup()
