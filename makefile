test:
	py.test -v  ./tests --cov=glycopeptidepy --cov-report html --cov-report term

retest:
	py.test -v ./tests --lf

update-cv-lists:
	python -m cogapp -r glycopeptidepy/io/cv/peff.py
	python -m autopep8 -i --max-line-length 80 glycopeptidepy/io/cv/peff.py

dev:
	python setup.py develop
