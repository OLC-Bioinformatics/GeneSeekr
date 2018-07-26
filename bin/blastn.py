#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import printtime
from methods.geneseekr import GeneSeekr, Parser, sequencenames
import multiprocessing


class BLASTn(object):

    def geneseekr(self):
        printtime('Performing {program} analyses on {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
        # Create the GeneSeekr object
        geneseekr = GeneSeekr()
        # Make blast databases (if necessary)
        printtime('Creating {at} blast databases as required'
                  .format(at=self.analysistype),
                  self.start)
        geneseekr.makeblastdb(self.combinedtargets)
        # Populate variables
        self.targetfolders, self.targetfiles, self.records = geneseekr.target_folders(self.metadata,
                                                                                      self.analysistype)
        printtime('Performing {program} analyses on {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
        self.metadata = geneseekr.run_blast(self.metadata,
                                            self.analysistype,
                                            self.program,
                                            self.outfmt,
                                            evalue=self.evalue,
                                            num_threads=self.cpus)
        # Parse the output depending on whether unique results are desired
        printtime('Parsing {program} results for {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
        if self.unique:
            # Run the unique blastn parsing module
            self.metadata = geneseekr.unique_parse_blast(self.metadata,
                                                         self.analysistype,
                                                         self.fieldnames,
                                                         self.cutoff,
                                                         self.program)
            # Filter the unique hits
            self.metadata = geneseekr.filter_unique(self.metadata,
                                                    self.analysistype)
        else:
            # Run the standard blastn parsing module
            self.metadata = geneseekr.parse_blastn(self.metadata,
                                                   self.analysistype,
                                                   self.fieldnames,
                                                   self.cutoff)
        # Create reports
        printtime('Creating {at} reports'.format(at=self.analysistype), self.start)
        if self.analysistype == 'resfinder':
            # ResFinder-specific report
            self.metadata = geneseekr.resfinder_reporter(self.metadata,
                                                         self.analysistype,
                                                         self.targetfolders,
                                                         self.reportpath,
                                                         self.align,
                                                         self.targetfiles,
                                                         self.records,
                                                         self.program)
        elif self.analysistype == 'virulence':
            # VirulenceFinder-specific report
            geneseekr.virulencefinder_reporter(self.metadata,
                                               self.analysistype,
                                               self.reportpath)
        else:
            # GeneSeekr-specific report
            self.metadata = geneseekr.reporter(self.metadata,
                                               self.analysistype,
                                               self.reportpath,
                                               self.align,
                                               self.targetfiles,
                                               self.records,
                                               self.program)
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata:
            delattr(sample[self.analysistype], "targetnames")
            delattr(sample[self.analysistype], "targets")
        printtime('{at} analyses complete'.format(at=self.analysistype), self.start)

    def __init__(self, args):
        self.cutoff = args.cutoff
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.numthreads if args.numthreads else multiprocessing.cpu_count() - 1
        self.align = args.align
        self.resfinder = args.resfinder
        self.virulencefinder = args.virulencefinder
        # Automatically set self.unique to true for ResFinder or VirulenceFinder analyses
        self.unique = True if self.resfinder or self.virulencefinder else args.unique
        self.analysistype = args.analysistype
        self.start = args.start
        self.evalue = args.evalue
        self.program = args.program
        # Run the Parser class from the geneseekr methods script to create lists of the database targets, and
        # combined targets, fasta sequences, and metadata objects.
        parse = Parser(self, args)
        parse.strainer()
        # Extract the variables from the object
        self.reportpath = parse.reportpath
        self.targets = parse.targets
        self.strains = parse.strains
        self.combinedtargets = parse.combinedtargets
        self.metadata = parse.metadata
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']
        self.outfmt = "'6 qseqid sseqid positive mismatch gaps " \
                      "evalue bitscore slen length qstart qend qseq sstart send sseq'"
        self.targetfolders = set()
        self.targetfiles = list()
        self.records = dict()
