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
install:
	cp $(RESULTS) $(patsubst %,../%,$(RESULTS))

clean:
	rm -f -r ./build
	rm -f *.so