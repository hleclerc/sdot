PYTHON = python

all: link_install

pip_upload:
	# sudo apt install python3-venv
	# ${PYTHON} -m pip install --upgrade --break-system-packages twine
	${PYTHON} -m build
	${PYTHON} -m twine upload --repository pypi --verbose dist/*

link_install:
	${PYTHON} -m pip install --break-system-packages -e .


.PHONY: all \
	pip_upload \
	link_install