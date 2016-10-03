test:
	nosetests --with-coverage --with-timer --cover-package=glycopeptidepy --cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id glycopeptidepy/test/

retest:
	nosetests --cover-package=glycopeptidepy --logging-level=DEBUG -v --with-id --failed glycopeptidepy/test/