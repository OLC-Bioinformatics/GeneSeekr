#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import printtime
from geneseekr.geneseekr import Fields, GeneSeekr


class tBLASTn(Fields):

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
            # Run the unique blast parsing module
            self.metadata = geneseekr.unique_parse_blast(self.metadata,
                                                         self.analysistype,
                                                         self.fieldnames,
                                                         self.cutoff,
                                                         self.program)
            # Filter the unique hits
            self.metadata = geneseekr.filter_unique(self.metadata,
                                                    self.analysistype)
        else:
            # Run the standard blast parsing module
            self.metadata = geneseekr.parse_blast(self.metadata,
                                                  self.analysistype,
                                                  self.fieldnames,
                                                  self.cutoff,
                                                  self.program)

        # Create dictionaries
        self.metadata = geneseekr.dict_initialise(self.metadata,
                                                  self.analysistype)
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
        args.program = 'tblastn'
        super().__init__(args)
