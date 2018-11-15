#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import combinetargets, GenObject, make_path, MetadataObject
from Bio.Sequencing.Applications import SamtoolsFaidxCommandline
from io import StringIO
from glob import glob
import logging
import os

__author__ = 'adamkoziol'


class Parser(object):

    def main(self):
        """
        Run the parsing methods
        """
        if not self.genus_specific:
            self.target_find()
            self.strainer()
        self.metadata_populate()

    def strainer(self):
        """
        Locate all the FASTA files in the supplied sequence path. Create basic metadata objects for
        each sample
        """
        assert os.path.isdir(self.sequencepath), 'Cannot locate sequence path as specified: {}' \
            .format(self.sequencepath)
        # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
        # begins with .fa
        self.strains = sorted(glob(os.path.join(self.sequencepath, '*.fa*'.format(self.sequencepath))))
        # Populate the metadata object. This object will be populated to mirror the objects created in the
        # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
        # of a pipeline
        assert self.strains, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                             'to ensure that your sequence path is correct'.format(self.sequencepath)
        for sample in self.strains:
            # Create the object
            metadata = MetadataObject()
            # Set the base file name of the sequence. Just remove the file extension
            filename = os.path.splitext(os.path.split(sample)[1])[0]
            # Set the .name attribute to be the file name
            metadata.name = filename
            # Create the .general attribute
            metadata.general = GenObject()
            # Set the .general.bestassembly file to be the name and path of the sequence file
            metadata.general.bestassemblyfile = sample
            # Append the metadata for each sample to the list of samples
            self.metadata.append(metadata)

    def target_find(self):
        """
        Locate all .tfa FASTA files in the supplied target path. If the combinedtargets.fasta file
        does not exist, run the combine targets method
        """
        self.targets = sorted(glob(os.path.join(self.targetpath, '*.tfa')))
        try:
            self.combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
        except IndexError:
            combinetargets(self.targets, self.targetpath)
            self.combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
        assert self.combinedtargets, 'Could not find any files with an extension starting with "fa" in {}. ' \
                                     'Please check to ensure that your target path is correct'.format(self.targetpath)

    def genus_targets(self, metadata):
        """
        Find all target files (files with .tfa extensions), and create a combined targets file if it already doesn't
        exist
        """
        if self.analysistype != 'GDCS':
            metadata[self.analysistype].targetpath = os.path.join(self.targetpath, metadata.general.referencegenus)
        else:
            metadata[self.analysistype].targetpath = os.path.join(self.targetpath,
                                                                  'GDCS',
                                                                  metadata.general.referencegenus)
        metadata[self.analysistype].targets = \
            sorted(glob(os.path.join(metadata[self.analysistype].targetpath, '*.tfa')))
        try:
            metadata[self.analysistype].combinedtargets = \
                glob(os.path.join(metadata[self.analysistype].targetpath, '*.fasta'))[0]
        except IndexError:
            try:
                combinetargets(self.targets, self.targetpath)
                metadata[self.analysistype].combinedtargets = \
                    glob(os.path.join(metadata[self.analysistype].targetpath, '*.fasta'))[0]
            except IndexError:
                metadata[self.analysistype].combinedtargets = 'NA'
        metadata[self.analysistype].targetnames = [os.path.splitext(os.path.basename(fasta))[0] for fasta in
                                                   metadata[self.analysistype].targets]

    def metadata_populate(self):
        """
        Populate the :analysistype GenObject
        """
        logging.info('Extracting sequence names from combined target file')
        # Extract all the sequence names from the combined targets file
        if not self.genus_specific:
            sequence_names = sequencenames(self.combinedtargets)
        else:
            sequence_names = list()
        for metadata in self.metadata:
            # Create and populate the :analysistype attribute
            setattr(metadata, self.analysistype, GenObject())
            if not self.genus_specific:
                metadata[self.analysistype].targets = self.targets
                metadata[self.analysistype].combinedtargets = self.combinedtargets
                metadata[self.analysistype].targetpath = self.targetpath
                metadata[self.analysistype].targetnames = sequence_names
            else:
                self.genus_targets(metadata)
            try:
                metadata[self.analysistype].reportdir = os.path.join(metadata.general.outputdirectory,
                                                                     self.analysistype)
            except (AttributeError, KeyError):
                metadata[self.analysistype].reportdir = self.reportpath

    def __init__(self, args, pipeline=False):
        self.analysistype = args.analysistype
        self.sequencepath = os.path.abspath(os.path.join(args.sequencepath))
        self.targetpath = os.path.abspath(os.path.join(args.targetpath))
        self.pipeline = pipeline
        if self.pipeline:
            if 'assembled' in self.targetpath or 'mlst' in self.targetpath.lower():
                self.targetpath = self.targetpath.rstrip('_assembled')
                # self.targetpath = self.targetpath.rstrip('_full')
                if 'rmlst' in self.targetpath:
                    self.targetpath = os.path.join(os.path.dirname(self.targetpath), 'rMLST')
                elif 'mlst' in self.targetpath:
                    self.targetpath = os.path.join(os.path.dirname(self.targetpath), 'MLST')
        assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}' \
            .format(self.targetpath)
        self.reportpath = os.path.abspath(os.path.join(args.reportpath))
        make_path(self.reportpath)
        assert os.path.isdir(self.reportpath), 'Cannot locate report path as specified: {}' \
            .format(self.reportpath)
        self.logfile = os.path.join(self.sequencepath, 'log.txt')
        try:
            self.metadata = args.metadata
        except AttributeError:
            self.metadata = list()
        self.strains = list()
        self.targets = list()
        self.combinedtargets = list()
        self.genus_specific = args.genus_specific


def sequencenames(contigsfile):
    """
    Takes a multifasta file and returns a list of sequence names
    :param contigsfile: multifasta of all sequences
    :return: list of all sequence names
    """
    sequences = list()
    fai_file = contigsfile + '.fai'
    if not os.path.isfile(fai_file):
        logging.info('Creating .fai file for {contigs}'.format(contigs=contigsfile))
        samindex = SamtoolsFaidxCommandline(reference=contigsfile)
        # Run the sam index command
        stdoutindex, stderrindex = map(StringIO, samindex(cwd=os.path.dirname(contigsfile)))
        logging.debug(stdoutindex.getvalue())
        logging.debug(stderrindex.getvalue())
    # Read in the sequence names
    with open(fai_file, 'r') as seq_file:
        for line in seq_file:
            allele = line.split('\t')[0]
            sequences.append(allele)
    return sequences


def objector(kw_dict, start):
    metadata = MetadataObject()
    for key, value in kw_dict.items():
        setattr(metadata, key, value)
    # Ensure that only a single analysis is specified
    analysis_count = 0
    # Set the analysis type based on the arguments provided
    if metadata.resfinder is True:
        metadata.analysistype = 'resfinder'
        analysis_count += 1
    elif metadata.virulencefinder is True:
        metadata.analysistype = 'virulence'
        analysis_count += 1
    elif metadata.mlst is True:
        metadata.analysistype = 'mlst'
        analysis_count += 1
    # Warn that only one type of analysis can be performed at a time
    elif analysis_count > 1:
        logging.warning('Cannot perform MLST, ResFinder and VirulenceFinder simultaneously. Please choose only one '
                        'of the -M -R and -V flags')
    # Default to GeneSeekr
    else:
        metadata.analysistype = 'geneseekr'
    # Add the start time variable to the object
    metadata.start = start
    return metadata, False


# noinspection PyProtectedMember
def modify_usage_error(subcommand):
    """
    Method to append the help menu to a modified usage error when a subcommand is specified, but options are missing
    """
    import click
    from click._compat import get_text_stderr
    from click.utils import echo

    def show(self, file=None):
        import sys
        if file is None:
            file = get_text_stderr()
        color = None
        if self.ctx is not None:
            color = self.ctx.color
        echo('Error: %s\n' % self.format_message(), file=file, color=color)
        # Set the sys.argv to be the first two arguments passed to the script if the subcommand was specified
        arg2 = sys.argv[1] if sys.argv[1] in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'] else str()
        sys.argv = [' '.join([sys.argv[0], arg2])] if arg2 else [sys.argv[0]]
        # Call the help
        subcommand(['--help'])

    click.exceptions.UsageError.show = show
