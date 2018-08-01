#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="geneseekr",
    version="0.0.1",
    packages=find_packages(),
    scripts=[
	'geneseekr/geneseekr',
	'methods/geneseekr.py',
	'bin/blastn.py',
	'bin/blastp.py',
	'bin/blastx.py',
	'bin/tblastn.py',
	'bin/tblastx.py'
	],
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/OLC-Bioinformatics/GeneSeekr",
    install_requires=[
	'biopython',
	'click',
	'numpy',	
	'olctools', 
	'xlsxwriter'
	]
)
