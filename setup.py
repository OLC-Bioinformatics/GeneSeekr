#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="geneseekr",
    version="0.2.0",
    packages=find_packages(),
    scripts=['GeneSeekr'],
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/OLC-Bioinformatics/GeneSeekr",
    install_requires=[
        'biopython',
        'click',
        'numpy',
        'olctools',
        'xlsxwriter',
        'pytest']
)
