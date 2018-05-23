# Makefile for numerov_c extension

RESULTS = $(wildcard *.so)

build: clean
build:
	python3 setup.py build_ext --inplace

install: build
install:
	cp $(RESULTS) $(patsubst %,../%,$(RESULTS))

clean:
	rm -f -r ./build
	rm -f *.so