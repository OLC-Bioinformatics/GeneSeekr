#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import MetadataObject
from MLSTsippr.mlst import GeneSippr as MLSTSippr
from sipprverse_reporter.reports import Reports
from geneseekr.geneseekr import GeneSeekr
from geneseekr.parser import Parser
import multiprocessing
from glob import glob
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

    def run_blast(self):
        """
        Perform BLAST analyses
        """
        logging.info('Performing {program} analyses on {at} targets'.format(program=self.program,
                                                                            at=self.analysistype))
        if 'mlst' in self.analysistype.lower():
            self.metadata = self.geneseekr.run_blast(self.metadata,
                                                     self.analysistype,
                                                     self.program,
                                                     self.outfmt,
                                                     evalue='1E-10',
                                                     num_threads=self.cpus,
                                                     perc_identity=self.cutoff)
        elif 'sixteens' in self.analysistype:
            self.metadata = self.geneseekr.run_blast(self.metadata,
                                                     self.analysistype,
                                                     self.program,
                                                     self.outfmt,
                                                     evalue='1E-20',
                                                     num_threads=self.cpus,
                                                     num_alignments=5000,
                                                     perc_identity=99)
        elif self.analysistype == 'GDCS':
            self.metadata = self.geneseekr.run_blast(self.metadata,
                                                     self.analysistype,
                                                     self.program,
                                                     self.outfmt,
                                                     evalue=self.evalue,
                                                     num_threads=self.cpus,
                                                     task='blastn')
        else:
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
        if 'sixteens' in self.analysistype:
            self.metadata = self.geneseekr.sixteens_parser(metadata=self.metadata,
                                                           analysistype=self.analysistype,
                                                           fieldnames=self.fieldnames,
                                                           cutoff=self.cutoff,
                                                           program=self.program)
        elif self.unique:
            # Run the unique blast parsing module
            self.metadata = self.geneseekr.unique_parse_blast(metadata=self.metadata,
                                                              analysistype=self.analysistype,
                                                              fieldnames=self.fieldnames,
                                                              cutoff=self.cutoff,
                                                              program=self.program)
            # Filter the unique hits
            self.metadata = self.geneseekr.filter_unique(metadata=self.metadata,
                                                         analysistype=self.analysistype)
        else:
            # Run the standard blast parsing module
            self.metadata = self.geneseekr.parse_blast(metadata=self.metadata,
                                                       analysistype=self.analysistype,
                                                       fieldnames=self.fieldnames,
                                                       cutoff=self.cutoff,
                                                       program=self.program)

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
            self.targetfolders, self.targetfiles, self.records = \
                self.geneseekr.target_folders(metadata=self.metadata,
                                              analysistype=self.analysistype)
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
        elif 'mlst' in self.analysistype.lower():
            # Adjust the analysis type to be consistently lowercase
            self.analysistype = self.analysistype.lower()
            # Create the necessary attributes for the MLST-typing method
            for sample in self.metadata:
                sample[self.analysistype].alleles = sorted(list(set(allele.split('_')[0]
                                                                    for allele in sample[self.analysistype]
                                                                    .targetnames)))
                # In order to work with the Enterobase cgMLST scheme that has underscores in the gene names (e.g.
                # AEJV01_03887, check for multiple underscores in the allele name, and act appropriately
                if len(sample[self.analysistype].alleles) > 53:
                    allele_set = set()
                    for allele in sample[self.analysistype].targetnames:
                        if len(allele.split('_')) == 3:
                            allele = '_'.join([allele.split('_')[0], allele.split('_')[1]])
                        else:
                            allele = allele.split('_')[0]
                        allele_set.add(allele)
                    sample[self.analysistype].alleles = sorted(list(allele_set))
                # Allele names attribute is apparently the same as the alleles attribute
                sample[self.analysistype].allelenames = sample[self.analysistype].alleles
                try:
                    sample[self.analysistype].profile = glob(os.path.join(sample[self.analysistype].targetpath,
                                                                          '*.txt'))[0]
                except IndexError:
                    sample[self.analysistype].profile = 'NA'
            # Create the typing object
            typing = MLST(args=self,
                          pipelinecommit='',
                          startingtime=self.start,
                          scriptpath='',
                          analysistype=self.analysistype,
                          cutoff=0.99,
                          pipeline=self.pipeline)
            # Perform typing, and create reports
            typing.reporter()
        elif 'sixteens' in self.analysistype:
            self.metadata = self.geneseekr.sixteens_reporter(metadata=self.metadata,
                                                             analysistype=self.analysistype,
                                                             reportpath=self.reportpath)
        elif self.analysistype == 'GDCS':
            self.metadata = self.geneseekr.gdcs_reporter(metadata=self.metadata,
                                                         analysistype=self.analysistype,
                                                         reportpath=self.reportpath)
        elif 'sero' in self.analysistype:
            self.metadata = self.geneseekr.sero_reporter(metadata=self.metadata,
                                                         analysistype=self.analysistype,
                                                         reportpath=self.reportpath)
        elif self.analysistype == 'genesippr':
            sippr_report = GenesipprReport(start=self.start,
                                           metadata=self.metadata,
                                           sequencepath=self.sequencepath,
                                           reportpath=self.reportpath,
                                           analysistype=self.analysistype)
            sippr_report.main()
        else:
            # GeneSeekr-specific report
            self.targetfolders, self.targetfiles, self.records = \
                self.geneseekr.target_folders(metadata=self.metadata,
                                              analysistype=self.analysistype)
            self.metadata = self.geneseekr.reporter(metadata=self.metadata,
                                                    analysistype=self.analysistype,
                                                    reportpath=self.reportpath,
                                                    align=self.align,
                                                    targetfiles=self.targetfiles,
                                                    records=self.records,
                                                    program=self.program)

    # noinspection PyNoneFunctionAssignment
    def clean_object(self):
        """
        Remove certain attributes from the object; they take up too much room on the .json report
        """
        self.metadata = self.geneseekr.clean_object(self.metadata,
                                                    self.analysistype)

    def __init__(self, args, analysistype='geneseekr', cutoff=70, program='blastn', genus_specific=False, unique=False,
                 evalue='1E-05', pipeline=True):
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
                self.analysistype = args.analysistype.lower()
            except AttributeError:
                self.analysistype = analysistype.lower()
                args.analysistype = analysistype.lower()
        elif analysistype == 'GDCS':
            self.analysistype = analysistype
        else:
            self.analysistype = analysistype.lower()
        try:
            self.resfinder = args.resfinder
        except AttributeError:
            self.resfinder = False
        try:
            self.virulencefinder = args.virulencefinder
        except AttributeError:
            self.virulencefinder = False
        try:
            self.typing = args.typing
        except AttributeError:
            self.typing = False
        # Automatically set self.unique to true for the appropriate analyses
        try:
            self.unique = True if unique or self.resfinder is True or self.virulencefinder is True or 'resfinder' in \
                                  self.analysistype or 'virulence' in self.analysistype or self.typing is True \
                                  or 'mlst' in self.analysistype.lower() or 'sixteens' in self.analysistype \
                                  else args.unique
        except AttributeError:
            self.unique = unique
        try:
            self.start = args.start
        except AttributeError:
            self.start = args.starttime
        try:
            self.evalue = args.evalue
        except AttributeError:
            self.evalue = evalue
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
        self.pipeline = pipeline
        try:
            self.metadata = args.runmetadata.samples
            parse = Parser(args=self,
                           pipeline=self.pipeline)
            if not self.genus_specific:
                parse.target_find()
            parse.metadata_populate()
        except (AttributeError, KeyError):
            # Run the Parser class from the GeneSeekr methods script to create lists of the database targets, and
            # combined targets, fasta sequences, and metadata objects.
            parse = Parser(args=self,
                           pipeline=self.pipeline)
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
        # Class variables required to use MLST
        self.path = self.sequencepath
        self.logfile = os.path.join(self.path, 'log')
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = self.metadata


class MLST(MLSTSippr):
    def reporter(self):
        analysistype = 'rmlst' if 'rmlst' in self.analysistype.lower() else 'mlst'
        # Populate self.plusdict in order to reuse parsing code from an assembly-based method
        for sample in self.runmetadata.samples:
            self.plusdict[sample.name] = dict()
            self.matchdict[sample.name] = dict()
            if sample.general.bestassemblyfile != 'NA':
                for gene in sample[analysistype].allelenames:
                    self.plusdict[sample.name][gene] = dict()
                    for allele, percentidentity in sample[analysistype].blastresults.items():
                        if gene in allele:
                            # Split the allele number from the gene name using the appropriate delimiter
                            if '_' in allele:
                                splitter = '_'
                            elif '-' in allele:
                                splitter = '-'
                            else:
                                splitter = ''
                            self.matchdict[sample.name].update({gene: allele.split(splitter)[-1]})
                            try:
                                self.plusdict[sample.name][gene][allele.split(splitter)[-1]][percentidentity] \
                                    = 10
                            except KeyError:
                                self.plusdict[sample.name][gene][allele.split(splitter)[-1]] = dict()
                                self.plusdict[sample.name][gene][allele.split(splitter)[-1]][percentidentity] \
                                    = 10
                    if gene not in self.matchdict[sample.name]:
                        self.matchdict[sample.name].update({gene: 'N'})
        self.profiler()
        self.sequencetyper()
        self.mlstreporter()


class GenesipprReport(object):

    def main(self):
        # Update the genesippr GenObject to include attributes necessary for the Genesippr reporter module
        for sample in self.runmetadata.samples:
            sample[self.analysistype].avgdepth = dict()
            sample[self.analysistype].standarddev = dict()
            if sample[self.analysistype].blastresults != 'NA':
                sample[self.analysistype].results = sample[self.analysistype].blastresults
            else:
                sample[self.analysistype].results = dict()
        report = Reports(inputobject=self)
        report.reporter()

    def __init__(self, start, metadata, sequencepath, reportpath, analysistype='genesippr'):
        self.starttime = start
        self.metadata = metadata
        self.path = sequencepath
        self.reportpath = reportpath
        self.logfile = os.path.join(self.path, 'log')
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = self.metadata
        self.analysistype = analysistype
