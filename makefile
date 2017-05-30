test:
	py.test -v  glycopeptidepy --cov=glycopeptidepy --cov-report html --cov-report term

retest:
	py.test -v glycopeptidepy --lf
