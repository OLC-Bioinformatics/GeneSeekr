#!/usr/bin/env python3
from geneseekr.geneseekr import GeneSeekr
from geneseekr.parser import Parser
import multiprocessing
import logging
import os

__author__ = 'adamkoziol'


class BLAST(object):

    def seekr(self):
        """
        Run the methods in the proper order
        """
        self.blast_db()
        self.run_blast()
        self.parse_results()
        self.create_reports()
        self.clean_object()
        logging.info('{at} analyses complete'.format(at=self.analysistype))

    def blast_db(self):
        """
        Make blast databases (if necessary)
        """
        logging.info('Creating {at} blast databases as required'.format(at=self.analysistype))
        for sample in self.metadata:
            self.geneseekr.makeblastdb(sample[self.analysistype].combinedtargets,
                                       self.program)
            self.targetfolders, self.targetfiles, self.records = \
                self.geneseekr.target_folders(self.metadata,
                                              self.analysistype)

    def run_blast(self):
        """
        Perform BLAST analyses
        """
        logging.info('Performing {program} analyses on {at} targets'.format(program=self.program,
                                                                            at=self.analysistype))
        self.metadata = self.geneseekr.run_blast(self.metadata,
                                                 self.analysistype,
                                                 self.program,
                                                 self.outfmt,
                                                 evalue=self.evalue,
                                                 num_threads=self.cpus)

    def parse_results(self):
        """
        Parse the output depending on whether unique results are desired
        """
        logging.info('Parsing {program} results for {at} targets'.format(program=self.program,
                                                                         at=self.analysistype))
        if self.unique:
            # Run the unique blast parsing module
            self.metadata = self.geneseekr.unique_parse_blast(self.metadata,
                                                              self.analysistype,
                                                              self.fieldnames,
                                                              self.cutoff,
                                                              self.program)
            # Filter the unique hits
            self.metadata = self.geneseekr.filter_unique(self.metadata,
                                                         self.analysistype)
        else:
            # Run the standard blast parsing module
            self.metadata = self.geneseekr.parse_blast(self.metadata,
                                                       self.analysistype,
                                                       self.fieldnames,
                                                       self.cutoff,
                                                       self.program)

    def create_reports(self):
        """
        Create reports
        """
        # Create dictionaries
        self.metadata = self.geneseekr.dict_initialise(self.metadata,
                                                       self.analysistype)
        # Create reports
        logging.info('Creating {at} reports'.format(at=self.analysistype))
        if 'resfinder' in self.analysistype:
            # ResFinder-specific report
            self.metadata = self.geneseekr.resfinder_reporter(metadata=self.metadata,
                                                              analysistype=self.analysistype,
                                                              reportpath=self.reportpath,
                                                              align=self.align,
                                                              targetfiles=self.targetfiles,
                                                              records=self.records,
                                                              program=self.program,
                                                              targetpath=self.targetpath)
        elif 'virulence' in self.analysistype:
            # VirulenceFinder-specific report
            self.geneseekr.virulencefinder_reporter(self.metadata,
                                                    self.analysistype,
                                                    self.reportpath)
        else:
            # GeneSeekr-specific report
            self.metadata = self.geneseekr.reporter(self.metadata,
                                                    self.analysistype,
                                                    self.reportpath,
                                                    self.align,
                                                    self.targetfiles,
                                                    self.records,
                                                    self.program)

    # noinspection PyNoneFunctionAssignment
    def clean_object(self):
        """
        Remove certain attributes from the object; they take up too much room on the .json report
        """
        self.metadata = self.geneseekr.clean_object(self.metadata,
                                                    self.analysistype)

    def __init__(self, args, analysistype='geneseekr', cutoff=70, program='blastn', genus_specific=False):
        try:
            args.program = args.program
        except AttributeError:
            args.program = program
        self.program = args.program
        try:
            self.cutoff = args.cutoff
        except AttributeError:
            self.cutoff = cutoff
        try:
            self.cpus = args.numthreads if args.numthreads else multiprocessing.cpu_count() - 1
        except AttributeError:
            self.cpus = args.cpus
        try:
            self.align = args.align
        except AttributeError:
            self.align = True
        if analysistype == 'geneseekr':
            try:
                self.analysistype = args.analysistype
            except AttributeError:
                self.analysistype = analysistype
                args.analysistype = analysistype
        else:
            self.analysistype = analysistype
        try:
            self.resfinder = args.resfinder
        except AttributeError:
            self.resfinder = False
        try:
            self.virulencefinder = args.virulencefinder
        except AttributeError:
            self.virulencefinder = False
        # Automatically set self.unique to true for ResFinder or VirulenceFinder analyses
        self.unique = True if self.resfinder or self.virulencefinder or 'resfinder' in self.analysistype \
            or self.analysistype == 'virulencefinder' else args.unique
        try:
            self.start = args.start
        except AttributeError:
            self.start = args.starttime
        try:
            self.evalue = args.evalue
        except AttributeError:
            self.evalue = '1E-05'
        try:
            self.sequencepath = args.sequencepath
        except AttributeError:
            self.sequencepath = str()
        try:
            self.targetpath = os.path.join(args.reffilepath, self.analysistype)
        except AttributeError:
            self.targetpath = args.targetpath
        self.reportpath = args.reportpath
        self.genus_specific = genus_specific
        try:
            self.metadata = args.runmetadata.samples
            parse = Parser(self)
            if not self.genus_specific:
                parse.target_find()
            parse.metadata_populate()
        except (AttributeError, KeyError):
            # Run the Parser class from the GeneSeekr methods script to create lists of the database targets, and
            # combined targets, fasta sequences, and metadata objects.
            parse = Parser(self)
            parse.main()
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
        # Create the GeneSeekr object
        self.geneseekr = GeneSeekr()
