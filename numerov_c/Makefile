# Makefile for numerov_c extension

RESULTS=$(wildcard *.so)
CANUSETHREE=$(shell python3 --version 2>/dev/null)
PY=python3

ifeq (CANUSETHREE,)
	PY=python
endif

build: clean
build:
	$(PY) setup.py build_ext --inplace

install: build
install: test
install:
	make cpfiles

cpfiles:
	cp $(RESULTS) $(patsubst %,../%,$(RESULTS))

test:
	$(PY) unit_tests.py

clean:
	rm -f -r ./build
	rm -f *.so