#!/usr/bin/env python
from ARMI_Lt import ARMISeekr
from Bio.Blast.Applications import NcbiblastnCommandline
__author__ = 'mike knowles'


class GeneSeekr(ARMISeekr):

    def _blast(self, (fasta, db)):
        self.yeah()
        blastn = NcbiblastnCommandline('/usr/local/bin/blastn',
                                       query=fasta,
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

