#!/usr/bin/env python3
from geneseekr.geneseekr import BLAST


class tBLASTx(BLAST):

    def __init__(self, args, analysistype='geneseekr', cutoff=70):
        super().__init__(args, analysistype, cutoff)