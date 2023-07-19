
PYPI_PASSWORD := $(shell cat pypi_pass.txt)

default:
	python setup.py install


push_to_pypi:
	rm -fr dist
	python3 -m build
	twine upload -r pypi dist/* --user yvesmartindestaillades --password $(PYPI_PASSWORD)
	rm -fr dist