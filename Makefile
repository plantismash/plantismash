PYTHON_FILES=$(shell find . -name '*.py')

all: unit coverage

tags: ${PYTHON_FILES}
	ctags --python-kinds=-i ${PYTHON_FILES}

coverage: ${PYTHON_FILES}
	rm -rf cover
	mkdir cover
	nosetests -v --with-coverage --cover-html --cover-package="antismash"

unit: ${PYTHON_FILES}
	nosetests -v

integration: ${PYTHON_FILES}
	nosetests -v -m "(?:^|[\b_\./-])[Ii]ntegration"

.PHONY: coverage unit integration
