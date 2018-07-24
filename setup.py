#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="GeneSeekr",
    version="0.0.1",
    packages=find_packages(),
    scripts=['methods/geneseekr',
	     'bin/blastn'],
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/OLC-Bioinformatics/GeneSeekr",
    install_requires=['olctools', 
		      'biopython',
		      'xlsxwriter
	]
)
