#!/usr/bin/env python
from ARMI_Lt import ARMISeekr
from Bio.Blast.Applications import NcbiblastnCommandline
__author__ = 'mike knowles'


class Jackson(ARMISeekr):

    def _blast(self, (fasta, db)):
        blastn = NcbiblastnCommandline(query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen'",
                                       perc_identity=self.cutoff,
                                       num_descriptions=10000,
                                       num_alignments=10000)
        stdout, stderr = blastn()
        if stdout != '':
            return [[fasta, aln[0],
                     (lambda x: "{0.2d}%%".format(x) if x <= 100.0 else '+')(abs(float(aln[1]) / float(aln[2])))]
                    for aln in [hsp.split('\t')
                    for hsp in stdout.rstrip().split("\n")]
                    if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]

