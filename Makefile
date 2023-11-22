
PYPI_PASSWORD := $(shell cat ~/.pypi_pass.txt)

default:
	python setup.py install

test:
	python -m pytest --envfile env

push_to_pypi:
	rm -fr dist
	python -m build
	twine upload -r pypi dist/* --user yvesmartindestaillades --password $(PYPI_PASSWORD) --skip-existing
	rm -fr dist