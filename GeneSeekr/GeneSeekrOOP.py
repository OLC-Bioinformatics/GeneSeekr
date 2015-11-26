#!/usr/bin/env python
from ARMI_Lt import ARMISeekr
from Bio.Blast.Applications import NcbiblastnCommandline
from argparse import ArgumentParser
import os
__author__ = 'mike knowles'


class GeneSeekr(ARMISeekr):

    def _blast(self, (fasta, db)):
        self.yeah()
        blastn = NcbiblastnCommandline(query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen'",
                                       perc_identity=self.cutoff,
                                       num_descriptions=10000,
                                       num_alignments=10000)
        stdout, stderr = blastn()
        if stdout != '':
            return [[fasta, aln[0].split('_')[0], int(aln[0].split('_')[1])]
                    for aln in [hsp.split('\t')
                                for hsp in stdout.rstrip().split("\n")]
                    if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]


parent = ArgumentParser(add_help=False)
parent.add_argument('--version', action='version', version='%(prog)s v0.5')
parent.add_argument('-i', '--input', required=True, help='Specify input fasta folder')
parent.add_argument('-o', '--output', default=os.getcwd(), help='Specify output folder for csv and json')
parent.add_argument('-c', '--cutoff', type=int, help='Threshold for maximum unique bacteria for a single allele')
parent.add_argument('-t', '--threads', type=int, default=12, help='Specify number of threads')
parent.add_argument('-m', '--marker', required=True, help='Specify gene file in FASTA format')
