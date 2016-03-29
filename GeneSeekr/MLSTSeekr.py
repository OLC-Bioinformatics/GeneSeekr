#!/usr/bin/env python
from ARMI_Lt import ARMISeekr
from Bio.Blast.Applications import NcbiblastnCommandline
from argparse import ArgumentParser
import os
import multiprocessing
__author__ = 'mike knowles'


class GeneSeekr(ARMISeekr):

    def _blast(self, (fasta, db)):
        blastn = NcbiblastnCommandline(query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen qacc'",
                                       perc_identity=self.cutoff,
                                       num_descriptions=10000,
                                       num_alignments=10000)
        stdout, stderr = blastn()
        if stdout != '':
            return [[[gene], [int(allele), qacc]]
                    for sseqid, nident, slen, qacc in [hsp.split('\t')
                    for hsp in stdout.rstrip().split("\n")]
                    for gene, allele in [sseqid.split('_')]
                    if abs(float(nident) / float(slen)) >= self.cutoff/100.0]


parent = ArgumentParser(add_help=False)
parent.add_argument('--version', action='version', version='%(prog)s v0.5')
parent.add_argument('input', nargs='?', default=os.path.relpath(os.getcwd()), help='Specify input fasta folder')
parent.add_argument('-o', '--output', default=os.getcwd(), help='Specify output folder for csv and json')
parent.add_argument('-c', '--cutoff', type=int, help='Threshold for maximum unique bacteria for a single allele')
parent.add_argument('-t', '--threads', type=int, default=multiprocessing.cpu_count(), help='Specify number of threads')
parent.add_argument('-m', '--marker', required=True, help='Specify gene file in FASTA format')
parent.add_argument('--evalue', default=1e-7, type=float, help='BLAST evalue to use (default 1e-7)')
