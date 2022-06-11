test:
	py.test -v  glycopeptidepy --cov=glycopeptidepy --cov-report html --cov-report term

retest:
	py.test -v glycopeptidepy --lf

update-cv-lists:
	python -m cogapp -r glycopeptidepy/io/cv/peff.py glycopeptidepy/structure/modification/data/psimod.py glycopeptidepy/structure/modification/data/uniprot.py
	python -m autopep8 -i --max-line-length 80 glycopeptidepy/io/cv/peff.py glycopeptidepy/structure/modification/data/psimod.py glycopeptidepy/structure/modification/data/uniprot.py

dev:
	python setup.py develop
