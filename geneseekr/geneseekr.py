#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, run_subprocess
from accessoryFunctions.resistance import ResistanceNotes
from sipprverse_reporter.reports import Reports
from biotools.bbtools import kwargs_to_string
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastxCommandline, NcbiblastpCommandline, \
    NcbitblastnCommandline, NcbitblastxCommandline
from Bio.Application import ApplicationError
from Bio.Data.CodonTable import TranslationError
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from click import progressbar
from csv import DictReader
from glob import glob
import xlsxwriter
import operator
import logging
import csv
import sys
import os
import re

__author__ = 'adamkoziol'


class GeneSeekr(object):

    @staticmethod
    def makeblastdb(fasta, program='blastn', returncmd=False, **kwargs):
        """
        Wrapper for makeblastdb. Assumes that makeblastdb is an executable in your $PATH
        Makes blast database files from targets as necessary
        :param fasta: Input FASTA-formatted file
        :param program: BLAST program used
        :param returncmd: Boolean for if the makeblastdb command should be returned
        :param kwargs: Dictionary of optional arguments
        :return: Stdout, Stderr, makeblastdb command (if requested)
        """
        # Convert the options dictionary to a string
        options = kwargs_to_string(kwargs)
        # Set the dbtype appropriately
        if program == 'blastn' or program == 'tblastn' or program == 'tblastx':
            dbtype = 'nucl'
        else:
            dbtype = 'prot'
        # Remove the file extension from the file name
        output = os.path.splitext(fasta)[0]
        cmd = 'makeblastdb -in {fasta} -parse_seqids -max_file_sz 2GB -dbtype {dbtype} -out {output}{options}' \
            .format(fasta=fasta,
                    dbtype=dbtype,
                    output=output,
                    options=options)
        # Check if database already exists
        if not os.path.isfile('{}.nhr'.format(output)):
            if output != 'NA':
                logging.debug(cmd)
            out, err = run_subprocess(cmd)
        else:
            out = str()
            err = str()
        if returncmd:
            return out, err, cmd
        else:
            return out, err

    @staticmethod
    def target_folders(metadata, analysistype):
        """
        Create a set of all database folders used in the analyses
        :param metadata: Metadata object
        :param analysistype: Name of analysis type
        :return: Lists of all target folders and files used in the analyses. Dictionary of SeqIO
        parsed sequences
        """
        targetfolders = set()
        targetfiles = list()
        records = dict()
        for sample in metadata:
            if sample[analysistype].combinedtargets != 'NA':
                targetfolders.add(sample[analysistype].targetpath)
        for targetdir in targetfolders:
            # List comprehension to remove any previously created database files from list
            targetfiles = glob(os.path.join(targetdir, '*.fasta'))
            for targetfile in targetfiles:
                # Read the sequences from the target file to a dictionary
                records[targetfile] = SeqIO.to_dict(SeqIO.parse(targetfile, 'fasta'))
        return targetfolders, targetfiles, records

    def run_blast(self, metadata, analysistype, program, outfmt, evalue='1E-5', num_threads=12, num_alignments=1000000,
                  perc_identity=70, task='blastn'):
        """
        Runs BLAST on all the samples in the metadata object
        :param metadata: Metadata object
        :param analysistype: Name of analysis type
        :param program: BLAST program to use for the alignment
        :param outfmt; Custom fields to include in BLAST output
        :param evalue: e-value cut-off for BLAST analyses
        :param num_threads: Number of threads to use for BLAST analyses
        :param num_alignments: Number of alignments to perform in BLAST analyses
        :param perc_identity: Percent identity cutoff
        :param task: For short sequences being analysed with BLASTn, allow for the blastn-short parameter to be
        specified,
        :return: Updated metadata object
        """
        with progressbar(metadata) as bar:
            for sample in bar:
                # Run the BioPython BLASTn module with the genome as query, fasta (target gene) as db.
                make_path(sample[analysistype].reportdir)
                # Set the name and path of the BLAST report as reportdir/samplename_blastprogram.csv
                sample[analysistype].report = os.path.join(
                    sample[analysistype].reportdir, '{name}_{program}.csv'.format(name=sample.name,
                                                                                  program=program))
                # Check the size of the report (if it exists). If it has size 0, something went wrong on a previous
                # iteration of the script. Delete the empty file in preparation for another try
                try:
                    size = os.path.getsize(sample[analysistype].report)
                    # If a report was created, but no results entered - program crashed, or no sequences passed
                    # thresholds, remove the report, and run the blast analyses again
                    if size == 0:
                        os.remove(sample[analysistype].report)
                except FileNotFoundError:
                    pass
                # Split the extension from the file path
                db = os.path.splitext(sample[analysistype].combinedtargets)[0]
                # Create the command line argument using the appropriate BioPython BLAST wrapper
                if program == 'blastn':
                    blast = self.blastn_commandline(sample=sample,
                                                    analysistype=analysistype,
                                                    db=db,
                                                    evalue=evalue,
                                                    num_alignments=num_alignments,
                                                    num_threads=num_threads,
                                                    outfmt=outfmt,
                                                    perc_identity=perc_identity,
                                                    task=task)
                elif program == 'blastp':
                    blast = self.blastp_commandline(sample=sample,
                                                    analysistype=analysistype,
                                                    db=db,
                                                    evalue=evalue,
                                                    num_alignments=num_alignments,
                                                    num_threads=num_threads,
                                                    outfmt=outfmt)
                elif program == 'blastx':
                    blast = self.blastx_commandline(sample=sample,
                                                    analysistype=analysistype,
                                                    db=db,
                                                    evalue=evalue,
                                                    num_alignments=num_alignments,
                                                    num_threads=num_threads,
                                                    outfmt=outfmt)
                elif program == 'tblastn':
                    blast = self.tblastn_commandline(sample=sample,
                                                     analysistype=analysistype,
                                                     db=db,
                                                     evalue=evalue,
                                                     num_alignments=num_alignments,
                                                     num_threads=num_threads,
                                                     outfmt=outfmt)
                elif program == 'tblastx':
                    blast = self.tblastx_commandline(sample=sample,
                                                     analysistype=analysistype,
                                                     db=db,
                                                     evalue=evalue,
                                                     num_alignments=num_alignments,
                                                     num_threads=num_threads,
                                                     outfmt=outfmt)
                else:
                    blast = str()
                assert blast, 'Something went wrong, the BLAST program you provided ({program}) isn\'t supported'\
                    .format(program=program)
                # Save the blast command in the metadata
                sample[analysistype].blastcommand = str(blast)
                # Only run blast if the report doesn't exist
                if not os.path.isfile(sample[analysistype].report):
                    try:
                        blast()
                    except ApplicationError:
                        try:
                            os.remove(sample[analysistype].report)
                        except (IOError, ApplicationError):
                            pass
        # Return the updated metadata object
        return metadata

    @staticmethod
    def blastn_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt, perc_identity,
                           task='blastn'):
        # BLAST command line call. Note the high number of default alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        blastn = NcbiblastnCommandline(query=sample.general.bestassemblyfile,
                                       db=db,
                                       evalue=evalue,
                                       task=task,
                                       num_alignments=num_alignments,
                                       num_threads=num_threads,
                                       outfmt=outfmt,
                                       perc_identity=perc_identity,
                                       out=sample[analysistype].report)
        return blastn

    @staticmethod
    def blastx_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt):
        blastx = NcbiblastxCommandline(query=sample.general.bestassemblyfile,
                                       db=db,
                                       evalue=evalue,
                                       num_alignments=num_alignments,
                                       num_threads=num_threads,
                                       outfmt=outfmt,
                                       out=sample[analysistype].report)
        return blastx

    @staticmethod
    def blastp_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt):
        blastp = NcbiblastpCommandline(query=sample.general.bestassemblyfile,
                                       db=db,
                                       evalue=evalue,
                                       num_alignments=num_alignments,
                                       num_threads=num_threads,
                                       outfmt=outfmt,
                                       out=sample[analysistype].report)
        return blastp

    @staticmethod
    def tblastn_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt):
        # BLAST command line call. Note the high number of default alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        tblastn = NcbitblastnCommandline(query=sample.general.bestassemblyfile,
                                         db=db,
                                         evalue=evalue,
                                         num_alignments=num_alignments,
                                         num_threads=num_threads,
                                         outfmt=outfmt,
                                         out=sample[analysistype].report)
        return tblastn

    @staticmethod
    def tblastx_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt):
        # BLAST command line call. Note the high number of default alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        tblastx = NcbitblastxCommandline(query=sample.general.bestassemblyfile,
                                         db=db,
                                         evalue=evalue,
                                         num_alignments=num_alignments,
                                         num_threads=num_threads,
                                         outfmt=outfmt,
                                         out=sample[analysistype].report)
        return tblastx

    @staticmethod
    def parse_blast(metadata, analysistype, fieldnames, cutoff, program):
        """
        Parse the blast results, and store necessary data in dictionaries in metadata object
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param fieldnames: List of column names in BLAST report
        :param cutoff: Percent identity threshold
        :param program: BLAST program used in the analyses
        :return: Updated metadata object
        """
        for sample in metadata:
            # Initialise a dictionary to store all the target sequences
            sample[analysistype].targetsequence = dict()
            try:
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
                resultdict = dict()
                # Go through each BLAST result
                for row in blastdict:
                    # Create the subject length variable - if the sequences are DNA (e.g. blastn), use the subject
                    # length as usual; if the sequences are protein (e.g. tblastx), use the subject length / 3
                    if program == 'blastn' or program == 'blastp' or program == 'blastx':
                        subject_length = float(row['subject_length'])

                    else:
                        subject_length = float(row['subject_length']) / 3
                    # Calculate the percent identity and extract the bitscore from the row
                    # Percent identity is the (length of the alignment - number of mismatches) / total subject length
                    percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                             subject_length * 100))
                    # Remove unwanted pipes added to the name
                    target = row['subject_id'].lstrip('gb|').rstrip('|') if '|' in row['subject_id'] else \
                        row['subject_id']
                    # If the percent identity is greater than the cutoff
                    if percentidentity >= cutoff:
                        # Update the dictionary with the target and percent identity
                        resultdict.update({target: percentidentity})
                        # Determine if the orientation of the sequence is reversed compared to the reference
                        if int(row['subject_end']) < int(row['subject_start']):
                            # Create a sequence object using Biopython
                            seq = Seq(row['query_sequence'], IUPAC.unambiguous_dna)
                            # Calculate the reverse complement of the sequence
                            querysequence = str(seq.reverse_complement())
                        # If the sequence is not reversed, use the sequence as it is in the output
                        else:
                            querysequence = row['query_sequence']
                        # Add the sequence in the correct orientation to the sample
                        sample[analysistype].targetsequence[target] = querysequence
                    # Add the percent identity to the object
                    sample[analysistype].blastresults = resultdict
                # Populate missing results with 'NA' values
                if len(resultdict) == 0:
                    sample[analysistype].blastresults = 'NA'
            except FileNotFoundError:
                sample[analysistype].blastresults = 'NA'
        return metadata

    @staticmethod
    def sixteens_parser(metadata, analysistype, fieldnames, cutoff, program):
        """
        Custom 16S parsing scheme - will determine the genus of all the BLAST hits (max 1000 returned from BLAST.
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param fieldnames: List of column names in BLAST report
        :param cutoff: Percent identity threshold
        :param program: BLAST program used in the analyses
        :return: Updated metadata object
        """
        dbrecords = dict()
        for sample in metadata:
            try:
                # Load the NCBI 16S reference database as a dictionary
                dbrecords = SeqIO.to_dict(SeqIO.parse(sample[analysistype].combinedtargets, 'fasta'))
                break
            except AttributeError:
                pass
        for sample in metadata:
            # Initialise a dictionary to store all the target sequences
            sample[analysistype].targetsequence = dict()
            # Initialise a dictionary to store the number of times a genus is the best hit
            sample[analysistype].frequency = dict()
            try:
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
                resultdict = dict()
                # Go through each BLAST result
                for row in blastdict:
                    # Create the subject length variable - if the sequences are DNA (e.g. blastn), use the subject
                    # length as usual; if the sequences are protein (e.g. tblastx), use the subject length / 3
                    if program == 'blastn' or program == 'blastp' or program == 'blastx':
                        subject_length = float(row['subject_length'])

                    else:
                        subject_length = float(row['subject_length']) / 3
                    # Calculate the percent identity and extract the bitscore from the row
                    # Percent identity is the (length of the alignment - number of mismatches) / total subject length
                    percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                             subject_length * 100))
                    target = row['subject_id']
                    # Extract the genus name. Use the subject id as a key in the dict of the reference db.
                    # It will return the record e.g. gi|1018196593|ref|NR_136472.1| Escherichia marmotae
                    # strain HT073016 16S ribosomal RNA, partial sequence
                    # This full description can be manipulated to extract the genus e.g. Escherichia
                    genus = dbrecords[target].description.split('|')[-1].split()[0]
                    # Increment the number of times this genus was found, or initialise the dictionary with this
                    # genus the first time it is seen
                    try:
                        sample[analysistype].frequency[genus] += 1
                    except KeyError:
                        sample[analysistype].frequency[genus] = 1
                    try:
                        resultdict[dbrecords[target].description] += 1
                    except KeyError:
                        resultdict[dbrecords[target].description] = 1
                    # Sort the dictionary based on the number of times a genus is seen
                    sample[analysistype].sortedgenera = sorted(sample[analysistype].frequency.items(),
                                                               key=operator.itemgetter(1), reverse=True)
                    try:
                        # Extract the top result, and set it as the genus of the sample
                        sample[analysistype].genus = sample[analysistype].sortedgenera[0][0]
                        # Previous code relies on having the closest refseq genus, so set this as above
                        # sample.general.closestrefseqgenus = sample[analysistype].genus
                        # sample.general.referencegenus = sample[analysistype].genus
                    except IndexError:
                        # Populate attributes with 'NA'
                        sample[analysistype].sortedgenera = 'NA'
                        sample[analysistype].genus = 'NA'
                        # sample.general.closestrefseqgenus = 'NA'
                        # sample.general.referencegenus = 'NA'
                    # If the percent identity is greater than the cutoff
                    if percentidentity >= cutoff:
                        # Update the dictionary with the target and percent identity
                        resultdict.update({target: percentidentity})
                        # Determine if the orientation of the sequence is reversed compared to the reference
                        if int(row['subject_end']) < int(row['subject_start']):
                            # Create a sequence object using Biopython
                            seq = Seq(row['query_sequence'], IUPAC.unambiguous_dna)
                            # Calculate the reverse complement of the sequence
                            querysequence = str(seq.reverse_complement())
                        # If the sequence is not reversed, use the sequence as it is in the output
                        else:
                            querysequence = row['query_sequence']
                        # Add the sequence in the correct orientation to the sample
                        sample[analysistype].targetsequence[target] = querysequence
                # Add the percent identity to the object
                sample[analysistype].blastresults = resultdict
                # Populate missing results with 'NA' values
                if len(resultdict) == 0:
                    sample[analysistype].blastresults = dict()
            except FileNotFoundError:
                sample[analysistype].blastresults = dict()
        return metadata

    @staticmethod
    def unique_parse_blast(metadata, analysistype, fieldnames, cutoff, program):
        """
        Find the best BLAST hit at a location
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param fieldnames: List of column names in BLAST report
        :param cutoff: Percent identity threshold
        :param program: BLAST program used in the analyses
        :return: Updated metadata object
        """
        for sample in metadata:
            # Initialise a dictionary to store all the target sequences
            sample[analysistype].targetsequence = dict()
            sample[analysistype].queryranges = dict()
            sample[analysistype].querypercent = dict()
            sample[analysistype].queryscore = dict()
            sample[analysistype].results = dict()
            try:
                # Encountering the following error: # _csv.Error: field larger than field limit (131072)
                # According to https://stackoverflow.com/a/15063941, increasing the field limit should fix the issue
                csv.field_size_limit(sys.maxsize)
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
                # Go through each BLAST result
                for row in blastdict:
                    # Create the subject length variable - if the sequences are DNA (e.g. blastn), use the subject
                    # length as usual; if the sequences are protein (e.g. tblastx), use the subject length / 3
                    if program == 'blastn' or program == 'blastp' or program == 'blastx':
                        subject_length = float(row['subject_length'])
                    else:
                        subject_length = float(row['subject_length']) / 3
                    # Calculate the percent identity
                    # Percent identity is the (length of the alignment - number of mismatches) / total subject length
                    percentidentity = float('{:0.2f}'.format((float(row['positives'])) / subject_length * 100))
                    target = row['subject_id'].lstrip('gb|').rstrip('|') if '|' in row['subject_id'] else \
                        row['subject_id']
                    contig = row['query_id']
                    high = max([int(row['query_start']), int(row['query_end'])])
                    low = min([int(row['query_start']), int(row['query_end'])])
                    score = row['bit_score']
                    # Create new entries in the blast results dictionaries with the calculated variables
                    row['percentidentity'] = percentidentity
                    row['low'] = low
                    row['high'] = high
                    row['alignment_fraction'] = float('{:0.2f}'.format(float(float(row['alignment_length']) /
                                                                             subject_length * 100)))
                    # If the percent identity is greater than the cutoff
                    if percentidentity >= cutoff:
                        try:
                            sample[analysistype].results[contig].append(row)
                            # Boolean to store whether the list needs to be updated
                            append = True
                            # Iterate through all the ranges. If the new range is different than any of the ranges
                            # seen before, append it. Otherwise, update the previous ranges with the longer range as
                            # necessary e.g. [2494, 3296] will be updated to [2493, 3296] with [2493, 3293], and
                            # [2494, 3296] will become [[2493, 3296], [3296, 4132]] with [3296, 4132]
                            for spot in sample[analysistype].queryranges[contig]:
                                # Update the low value if the new low value is slightly lower than before
                                if 1 <= (spot[0] - low) <= 100:
                                    # Update the low value
                                    spot[0] = low
                                    # It is not necessary to append
                                    append = False
                                # Update the previous high value if the new high value is slightly higher than before
                                elif 1 <= (high - spot[1]) <= 100:
                                    # Update the high value in the list
                                    spot[1] = high
                                    # It is not necessary to append
                                    append = False
                                # Do not append if the new low is slightly larger than before
                                elif 1 <= (low - spot[0]) <= 100:
                                    append = False
                                # Do not append if the new high is slightly smaller than before
                                elif 1 <= (spot[1] - high) <= 100:
                                    append = False
                                # Do not append if the high and low are the same as the previously recorded values
                                elif low == spot[0] and high == spot[1]:
                                    append = False
                            # If the result appears to be in a new location, add the data to the object
                            if append:
                                sample[analysistype].queryranges[contig].append([low, high])
                                sample[analysistype].querypercent[contig] = percentidentity
                                sample[analysistype].queryscore[contig] = score
                        # Initialise and populate the dictionary for each contig
                        except KeyError:
                            sample[analysistype].queryranges[contig] = list()
                            sample[analysistype].queryranges[contig].append([low, high])
                            sample[analysistype].querypercent[contig] = percentidentity
                            sample[analysistype].queryscore[contig] = score
                            sample[analysistype].results[contig] = list()
                            sample[analysistype].results[contig].append(row)
                            sample[analysistype].targetsequence[target] = dict()
                        # Determine if the query sequence is in a different frame than the subject, and correct
                        # by setting the query sequence to be the reverse complement
                        if int(row['subject_end']) < int(row['subject_start']):
                            # Create a sequence object using Biopython
                            seq = Seq(row['query_sequence'], IUPAC.unambiguous_dna)
                            # Calculate the reverse complement of the sequence
                            querysequence = str(seq.reverse_complement())
                        # If the sequence is not reversed, use the sequence as it is in the output
                        else:
                            querysequence = row['query_sequence']
                        # Add the sequence in the correct orientation to the sample
                        sample[analysistype].targetsequence[target] = querysequence
            except FileNotFoundError:
                pass
        # Return the updated metadata object
        return metadata

    @staticmethod
    def filter_unique(metadata, analysistype):
        """
        Filters multiple BLAST hits in a common region of the genome. Leaves only the best hit
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :return: Updated metaata object
        """
        for sample in metadata:
            # Initialise variables
            sample[analysistype].blastresults = dict()
            sample[analysistype].blastlist = list()
            resultdict = dict()
            rowdict = dict()
            try:
                # Iterate through all the contigs, which had BLAST hits
                for contig in sample[analysistype].queryranges:
                    # Find all the locations in each contig that correspond to the BLAST hits
                    for location in sample[analysistype].queryranges[contig]:
                        # Extract the BLAST result dictionary for the contig
                        for row in sample[analysistype].results[contig]:
                            # Initialise variable to reduce the number of times row['value'] needs to be typed
                            contig = row['query_id']
                            high = row['high']
                            low = row['low']
                            percentidentity = row['percentidentity']
                            # Join the two ranges in the location list with a comma
                            locstr = ','.join([str(x) for x in location])
                            # Create a set of the location of all the base pairs between the low and high (-1) e.g.
                            # [6, 10] would give 6, 7, 8, 9, but NOT 10. This turns out to be useful, as there are
                            # genes located back-to-back in the genome e.g. strB and strA, with locations of 2557,3393
                            # and 3393,4196, respectively. By not including 3393 in the strB calculations, I don't
                            # have to worry about this single bp overlap
                            loc = set(range(low, high))
                            # Use a set intersection to determine whether the current result overlaps with location
                            # This will allow all the hits to be grouped together based on their location
                            if loc.intersection(set(range(location[0], location[1]))):
                                # Populate the grouped hits for each location
                                try:
                                    resultdict[contig][locstr].append(percentidentity)
                                    rowdict[contig][locstr].append(row)
                                # Initialise and populate the lists of the nested dictionary
                                except KeyError:
                                    try:
                                        resultdict[contig][locstr] = list()
                                        resultdict[contig][locstr].append(percentidentity)
                                        rowdict[contig][locstr] = list()
                                        rowdict[contig][locstr].append(row)
                                    # As this is a nested dictionary, it needs to be initialised here
                                    except KeyError:
                                        resultdict[contig] = dict()
                                        resultdict[contig][locstr] = list()
                                        resultdict[contig][locstr].append(percentidentity)
                                        rowdict[contig] = dict()
                                        rowdict[contig][locstr] = list()
                                        rowdict[contig][locstr].append(row)
            except KeyError:
                pass
            # Dictionary of results
            results = dict()
            # Find the best hit for each location based on percent identity
            for contig in resultdict:
                # Do not allow the same gene to be added to the dictionary more than once
                genes = list()
                for location in resultdict[contig]:
                    # Initialise a variable to determine whether there is already a best hit found for the location
                    multiple = False
                    # Iterate through the BLAST results to find the best hit
                    for row in rowdict[contig][location]:
                        # Add the best hit to the .blastresults attribute of the object
                        if row['percentidentity'] == max(resultdict[contig][location]) and not multiple \
                                and row['subject_id'] not in genes:
                            # Update the list with the blast results
                            sample[analysistype].blastlist.append(row)
                            results.update({row['subject_id']: row['percentidentity']})
                            genes.append(row['subject_id'])
                            multiple = True
            # Add the dictionary of results to the metadata object
            sample[analysistype].blastresults = results
        # Return the updated metadata object
        return metadata

    @staticmethod
    def dict_initialise(metadata, analysistype):
        """
        Initialise dictionaries for storing DNA and amino acid sequences
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :return: Updated metadata
        """
        for sample in metadata:
            sample[analysistype].dnaseq = dict()
            sample[analysistype].protseq = dict()
            sample[analysistype].ntindex = dict()
            sample[analysistype].aaindex = dict()
            sample[analysistype].ntalign = dict()
            sample[analysistype].aaalign = dict()
            sample[analysistype].aaidentity = dict()
        return metadata

    def reporter(self, metadata, analysistype, reportpath, align, targetfiles, records, program):
        """
        Custom reports for standard GeneSeekr analyses.
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param reportpath: Path of folder in which report is to be created
        :param align: Boolean of whether alignments between query and subject sequences are desired
        :param targetfiles: List of all files used in the analyses
        :param records: Dictionary of SeqIO parsed sequence records
        :param program: BLAST program used to perform analyses
        :return: Updated metadata object
        """
        # Also make a CSV file with different formatting for portal parsing purposes
        # Format as: Strain,Gene1,Gene2
        #            ID,PercentID,PercentID for all strains input - have a zero when gene wasn't found.
        csv_output = os.path.join(reportpath, '{at}_{program}.csv'.format(at=analysistype,
                                                                          program=program))
        targets = list()
        for record in records:
            for item in records[record]:
                targets.append(item)
        with open(csv_output, 'w') as outfile:
            outfile.write('Strain')
            for target in targets:
                outfile.write(',{}'.format(target))
            outfile.write('\n')
            for sample in metadata:
                outfile.write('{}'.format(sample.name))
                for target in targets:
                    if target in sample[analysistype].blastresults:
                        outfile.write(',{}'.format(sample[analysistype].blastresults[target]))
                    else:
                        outfile.write(',0')
                outfile.write('\n')
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook(os.path.join(reportpath, '{at}_{program}.xlsx'
                                                    .format(at=analysistype,
                                                            program=program)))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        # Initialise a list of all the headers with 'Strain'
        headers = ['Strain']
        # A set to store whether which genes (if any) require additional columns in the report due to
        # printing of alignments
        align_set = set()
        # Determine which genes require alignments
        for sample in metadata:
            if sample[analysistype].targetnames != 'NA':
                if sample[analysistype].blastresults != 'NA':
                    for target in sorted(sample[analysistype].targetnames):
                        try:
                            # Only if the alignment option is selected, for inexact results, add alignments
                            if align and sample[analysistype].blastresults[target] != 100.00:
                                # Add the target to the set of targets requiring alignments
                                align_set.add(target)
                        except KeyError:
                            pass
        # Create the headers as required for targets with alignments
        for sample in metadata:
            if sample[analysistype].targetnames != 'NA':
                if sample[analysistype].blastresults != 'NA':
                    for target in sorted(sample[analysistype].targetnames):
                        # Add the name of the gene to the header
                        headers.append(target)
                        if target in align_set:
                            if program == 'blastn':
                                # Add the appropriate headers
                                headers.extend(['{target}_FASTA'.format(target=target),
                                                '{target}_aa_Alignment'.format(target=target),
                                                '{target}_aa_SNP_location'.format(target=target),
                                                '{target}_nt_Alignment'.format(target=target),
                                                '{target}_nt_SNP_location'.format(target=target)
                                                ])
                            else:
                                headers.extend(['{target}_FASTA'.format(target=target),
                                                '{target}_aa_Alignment'.format(target=target),
                                                '{target}_aa_SNP_location'.format(target=target),
                                                ])
                    # Only need to iterate through this once
                    break
        # Set the column to zero
        col = 0
        # Write the header to the spreadsheet
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set the column width based on the longest header
            try:
                columnwidth[col] = len(header) if len(header) > columnwidth[col] else columnwidth[col]
            except KeyError:
                columnwidth[col] = len(header)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
        for sample in metadata:
            # Initialise a list to store all the data for each strain
            data = list()
            if sample[analysistype].targetnames != 'NA':
                # Append the sample name to the data list only if the script could find targets
                data.append(sample.name)
                if sample[analysistype].blastresults != 'NA':
                    for target in sorted(sample[analysistype].targetnames):
                        try:
                            # Only if the alignment option is selected, for inexact results, add alignments
                            if align and sample[analysistype].blastresults[target] != 100.00:
                                # Align the protein (and nucleotide) sequences to the reference
                                sample = self.alignprotein(sample, analysistype, target, targetfiles, records, program)
                                # Create a FASTA-formatted sequence output of the query sequence
                                if program == 'blastn':
                                    record = SeqRecord(sample[analysistype].dnaseq[target],
                                                       id='{}_{}'.format(sample.name, target),
                                                       description='')
                                else:
                                    record = SeqRecord(sample[analysistype].protseq[target],
                                                       id='{}_{}'.format(sample.name, target),
                                                       description='')

                                # Add the alignment, and the location of mismatches for both nucleotide and amino
                                # acid sequences
                                if program == 'blastn':
                                    data.extend([str(sample[analysistype].blastresults[target]),
                                                 record.format('fasta'),
                                                 sample[analysistype].aaalign[target],
                                                 sample[analysistype].aaindex[target],
                                                 sample[analysistype].ntalign[target],
                                                 sample[analysistype].ntindex[target]
                                                 ])
                                else:
                                    data.extend([str(sample[analysistype].blastresults[target]),
                                                 record.format('fasta'),
                                                 sample[analysistype].aaalign[target],
                                                 sample[analysistype].aaindex[target],
                                                 ])
                            elif align and sample[analysistype].blastresults[target] == 100.00:
                                if target in align_set:
                                    if program == 'blastn':
                                        data.extend([str(sample[analysistype].blastresults[target]), '-', '-', '-',
                                                     '-', '-'])
                                    else:
                                        data.extend([str(sample[analysistype].blastresults[target]), '-', '-', '-'])
                                else:
                                    data.append(str(sample[analysistype].blastresults[target]))
                            # elif not align and sample[analysistype].blastresults[target] == 100.00:
                            #     data.append(str(sample[analysistype].blastresults[target]))
                            elif not align and float(sample[analysistype].blastresults[target]) > 95.00:
                                data.append(str(sample[analysistype].blastresults[target]))
                            else:
                                if target in align_set:
                                    if program == 'blastn':
                                        data.extend(['-', '-', '-', '-', '-', '-'])
                                    else:
                                        data.extend(['-', '-', '-', '-'])
                                else:
                                    data.append('-')
                        # If there are no blast results for the target, add a '-'
                        except (KeyError, TypeError):
                            if target in align_set:
                                if program == 'blastn':
                                    data.extend(['-', '-', '-', '-', '-', '-'])
                                else:
                                    data.extend(['-', '-', '-', '-'])
                            else:
                                data.append('-')
                # If there are no blast results at all, add a '-'
                else:
                    data.extend(['-'] * (len(headers) - 1))

            # Increment the row and reset the column to zero in preparation of writing results
            row += 1
            col = 0
            # List of the number of lines for each result
            totallines = list()
            # Write out the data to the spreadsheet
            for results in data:
                worksheet.write(row, col, results, courier)
                try:
                    # Counting the length of multi-line strings yields columns that are far too wide, only count
                    # the length of the string up to the first line break
                    alignmentcorrect = len(results.split('\n')[0])
                    # Count the number of lines for the data
                    lines = results.count('\n') if results.count('\n') >= 1 else 1
                    # Add the number of lines to the list
                    totallines.append(lines)
                # If there are no newline characters, set the width to the length of the string
                except AttributeError:
                    alignmentcorrect = len(results)
                    lines = 1
                    # Add the number of lines to the list
                    totallines.append(lines)
                # Increase the width of the current column, if necessary
                try:
                    columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else columnwidth[col]
                except KeyError:
                    columnwidth[col] = alignmentcorrect
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
            # Set the width of the row to be the number of lines (number of newline characters) * 12
            if len(totallines) != 0:
                worksheet.set_row(row, max(totallines) * 15)
            else:
                worksheet.set_row(row, 1)
        # Close the workbook
        workbook.close()
        # Return the updated metadata object
        return metadata

    def resfinder_reporter(self, metadata, analysistype, reportpath, align, targetfiles, records,
                           program, targetpath):
        """
        Custom reports for ResFinder analyses. These reports link the gene(s) found to their resistance phenotypes
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param reportpath: Path of folder in which report is to be created
        :param align: Boolean of whether alignments between query and subject sequences are desired
        :param targetfiles: List of all files used in the analyses
        :param records: Dictionary of SeqIO parsed sequence records
        :param program: BLAST program used in the analyses
        :param targetpath: Name and path of the folder containing the targets
        :return: Updated metadata object
        """
        # Since the resfinder database is used for both sipping and assembled analyses, but the analysis type is
        # different, strip off the _assembled, so the targets are set correctly
        targetpath = targetpath if analysistype != 'resfinder_assembled' else targetpath.rstrip('_assembled')
        resistance_classes = ResistanceNotes.classes(targetpath)
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook(os.path.join(reportpath, '{at}_{program}.xlsx'
                                                    .format(at=analysistype,
                                                            program=program)))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        extended = False
        percentage = 'PercentIdentity' if program == 'blastn' else 'PercentPositive'
        headers = ['Strain', 'Gene', 'Allele', 'Resistance', percentage, 'PercentCovered', 'Contig', 'Location']
        # Add the appropriate string to the headers based on whether the BLAST outputs are DNA/amino acids
        headers.append('nt_sequence') if program == 'blastn' else headers.append('aa_sequence')
        for sample in metadata:
            sample[analysistype].sampledata = list()
            sample[analysistype].pipelineresults = dict()
            # Process the sample only if the script could find targets
            if sample[analysistype].blastlist != 'NA' and sample[analysistype].blastlist:
                for result in sample[analysistype].blastlist:
                    # Set the name to avoid writing out the dictionary[key] multiple times
                    name = result['subject_id']
                    try:
                        # Extract the necessary variables from the gene name string
                        gname, genename, accession, allele = ResistanceNotes.gene_name(name)
                    except ValueError:
                        genename = name
                        allele = str()
                    # Initialise a list to store all the data for each strain
                    data = list()
                    # Determine resistance phenotype of the gene
                    resistance = ResistanceNotes.resistance(name, resistance_classes)
                    # Append the necessary values to the data list
                    data.append(genename)
                    data.append(allele)
                    data.append(resistance)
                    percentid = result['percentidentity']
                    data.append(percentid)
                    data.append(result['alignment_fraction'])
                    data.append(result['query_id'])
                    data.append('...'.join([str(result['low']), str(result['high'])]))
                    # Populate the .pipelineresults attribute for compatibility with the assembly pipeline
                    if percentid > 80:
                        try:
                            if genename not in sample[analysistype].pipelineresults[resistance]:
                                sample[analysistype].pipelineresults[resistance]\
                                    .append('{rgene} ({pid}%)'.format(rgene=genename,
                                                                      pid=percentid))
                        except KeyError:
                            sample[analysistype].pipelineresults[resistance] = list()
                            if genename not in sample[analysistype].pipelineresults[resistance]:
                                sample[analysistype].pipelineresults[resistance]\
                                    .append('{rgene} ({pid}%)'.format(rgene=genename,
                                                                      pid=percentid))
                    try:
                        # Only if the alignment option is selected, for inexact results, add alignments
                        if align and percentid != 100.00:
                            # Align the protein (and nucleotide) sequences to the reference
                            sample = self.alignprotein(sample, analysistype, name, targetfiles, records, program)
                            if not extended:
                                if program == 'blastn':
                                    # Add the appropriate headers
                                    headers.extend(['aa_Identity',
                                                    'aa_Alignment',
                                                    'aa_SNP_location',
                                                    'nt_Alignment',
                                                    'nt_SNP_location'
                                                    ])
                                else:
                                    headers.extend(['aa_Identity',
                                                    'aa_Alignment',
                                                    'aa_SNP_location',
                                                    ])
                                extended = True
                            # Create a FASTA-formatted sequence output of the query sequence
                            if program == 'blastn':
                                record = SeqRecord(sample[analysistype].dnaseq[name],
                                                   id='{}_{}'.format(sample.name, name),
                                                   description='')
                            else:
                                record = SeqRecord(sample[analysistype].protseq[name],
                                                   id='{}_{}'.format(sample.name, name),
                                                   description='')

                            # Add the alignment, and the location of mismatches for both nucleotide and amino
                            # acid sequences
                            if program == 'blastn':
                                data.extend([record.format('fasta'),
                                             sample[analysistype].aaidentity[name],
                                             sample[analysistype].aaalign[name],
                                             sample[analysistype].aaindex[name],
                                             sample[analysistype].ntalign[name],
                                             sample[analysistype].ntindex[name]
                                             ])
                            else:
                                data.extend([record.format('fasta'),
                                             sample[analysistype].aaidentity[name],
                                             sample[analysistype].aaalign[name],
                                             sample[analysistype].aaindex[name],
                                             ])
                        else:
                            if program == 'blastn':
                                record = SeqRecord(Seq(result['query_sequence'], IUPAC.unambiguous_dna),
                                                   id='{}_{}'.format(sample.name, name),
                                                   description='')
                            else:
                                record = SeqRecord(Seq(result['query_sequence'], IUPAC.protein),
                                                   id='{}_{}'.format(sample.name, name),
                                                   description='')
                            data.append(record.format('fasta'))
                            if align:
                                # Add '-'s for the empty results, as there are no alignments for exact matches
                                data.extend(['-', '-', '-', '-', '-'])
                    # If there are no blast results for the target, add a '-'
                    except (KeyError, TypeError):
                        data.append('-')
                    sample[analysistype].sampledata.append(data)
        if 'nt_sequence' not in headers and program == 'blastn':
            headers.append('nt_sequence')
            # Write the header to the spreadsheet
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set the column width based on the longest header
            try:
                columnwidth[col] = len(header) if len(header) > columnwidth[col] else columnwidth[
                    col]
            except KeyError:
                columnwidth[col] = len(header)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
            # Increment the row and reset the column to zero in preparation of writing results
        row += 1
        col = 0
        # Write out the data to the spreadsheet
        for sample in metadata:
            if not sample[analysistype].sampledata:
                worksheet.write(row, col, sample.name, courier)
                # Increment the row and reset the column to zero in preparation of writing results
                row += 1
                col = 0
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row)
                worksheet.set_column(col, col, columnwidth[col])
            for data in sample[analysistype].sampledata:
                columnwidth[col] = len(sample.name) + 2
                worksheet.set_column(col, col, columnwidth[col])
                worksheet.write(row, col, sample.name, courier)
                col += 1
                # List of the number of lines for each result
                totallines = list()
                for results in data:
                    #
                    worksheet.write(row, col, results, courier)
                    try:
                        # Counting the length of multi-line strings yields columns that are far too wide, only count
                        # the length of the string up to the first line break
                        alignmentcorrect = len(str(results).split('\n')[1])
                        # Count the number of lines for the data
                        lines = results.count('\n') if results.count('\n') >= 1 else 1
                        # Add the number of lines to the list
                        totallines.append(lines)
                    except IndexError:
                        try:
                            # Counting the length of multi-line strings yields columns that are far too wide, only count
                            # the length of the string up to the first line break
                            alignmentcorrect = len(str(results).split('\n')[0])
                            # Count the number of lines for the data
                            lines = results.count('\n') if results.count('\n') >= 1 else 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                        # If there are no newline characters, set the width to the length of the string
                        except AttributeError:
                            alignmentcorrect = len(str(results))
                            lines = 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                    # Increase the width of the current column, if necessary
                    try:
                        columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else \
                            columnwidth[col]
                    except KeyError:
                        columnwidth[col] = alignmentcorrect
                    worksheet.set_column(col, col, columnwidth[col])
                    col += 1
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row, max(totallines) * 11)
                # Increase the row counter for the next strain's data
                row += 1
                col = 0
        # Close the workbook
        workbook.close()
        # Return the updated metadata object
        return metadata

    @staticmethod
    def virulencefinder_reporter(metadata, analysistype, reportpath):
        """
        Custom reports for VirulenceFinder analyses. These reports link the gene(s) found to their virulence phenotypes
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param reportpath: Path of folder in which report is to be created
        """
        with open(os.path.join(reportpath, 'virulence.csv'), 'w') as report:
            header = 'Strain,Gene,PercentIdentity,PercentCovered,Contig,Location,Sequence\n'
            data = ''
            for sample in metadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[analysistype].blastlist:
                        data += '{},'.format(sample.name)
                        multiple = False
                        for result in sample[analysistype].blastlist:
                            if analysistype == 'virulence':
                                gene = result['subject_id'].split(':')[0]
                            else:
                                gene = result['subject_id']
                            if multiple:
                                data += ','
                            data += '{},{},{},{},{}..{},{}\n' \
                                .format(gene, result['percentidentity'], result['alignment_fraction'],
                                        result['query_id'], result['low'], result['high'], result['query_sequence'])
                            # data += '\n'
                            multiple = True
                    else:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(header)
            report.write(data)

    @staticmethod
    def sixteens_reporter(metadata, analysistype, reportpath):
        """
        Custom reports for VirulenceFinder analyses. These reports link the gene(s) found to their virulence phenotypes
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param reportpath: Path of folder in which report is to be created
        """
        # Find the best 16S match
        for sample in metadata:
            if sample.general.bestassemblyfile != 'NA':
                if sample[analysistype].blastresults:
                    # Sort the dictionary based on the percent identity - set the highest results as the best hit
                    sample[analysistype].besthit = sorted(sample[analysistype].blastresults.items(),
                                                          key=operator.itemgetter(1), reverse=True)[0]
                else:
                    sample[analysistype].besthit = str()
            else:
                sample[analysistype].besthit = str()
        # Create the report
        with open(os.path.join(reportpath, 'sixteens.csv'), 'w') as report:
            header = 'Strain,Gene,PercentIdentity,Genus\n'
            data = ''
            for sample in metadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[analysistype].besthit:
                        sample[analysistype].sixteens_match = sample[analysistype].besthit[0].replace(',', '')
                        data += '{},'.format(sample.name)
                        multiple = False
                        if multiple:
                            data += ','
                        data += '{gene},{pi},{genus}\n' \
                            .format(gene=sample[analysistype].sixteens_match,
                                    pi=sample[analysistype].besthit[1],
                                    genus=sample.general.closestrefseqgenus)
                    else:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(header)
            report.write(data)
        # Return the updated metadata object
        return metadata

    @staticmethod
    def gdcs_reporter(metadata, analysistype, reportpath):
        """
        Creates a report of the GDCS results
        :param metadata: Metadata object
        :param analysistype: The variable to use when accessing attributes in the metadata object
        :param reportpath: Path of folder in which report is to be created
        """
        logging.info('Creating {} report'.format(analysistype))
        # Initialise list to store all the GDCS genes, and genera in the analysis
        gdcs = list()
        genera = list()
        for sample in metadata:
            sample[analysistype].faidict = dict()
            if sample.general.bestassemblyfile != 'NA':
                if os.path.isdir(sample[analysistype].targetpath):
                    # Update the fai dict with all the genes in the analysis, rather than just those with baited hits
                    Reports.gdcs_fai(sample)
                    sample[analysistype].createreport = True
                    # Determine which genera are present in the analysis
                    if sample.general.closestrefseqgenus not in genera:
                        genera.append(sample.general.closestrefseqgenus)
                    try:
                        # Add all the GDCS genes to the list
                        for gene in sorted(sample[analysistype].faidict):
                            if gene not in gdcs:
                                gdcs.append(gene)
                    except AttributeError:
                        sample[analysistype].createreport = False
                else:
                    sample[analysistype].createreport = False
            else:
                sample[analysistype].createreport = False
                sample.general.incomplete = True
        header = 'Strain,Genus,Matches,Pass/Fail,{},\n'.format(','.join(gdcs))
        data = str()
        with open(os.path.join(reportpath, '{}.csv'.format(analysistype)), 'w') as report:
            # Sort the samples in the report based on the closest refseq genus e.g. all samples with the same genus
            # will be grouped together in the report
            for genus in genera:
                for sample in metadata:
                    if sample.general.closestrefseqgenus == genus:
                        if sample[analysistype].createreport:
                            sample[analysistype].totaldepth = list()
                            # Add the sample to the report if it matches the current genus
                            # if genus == sample.general.closestrefseqgenus:
                            data += '{},{},'.format(sample.name, genus)
                            # Initialise a variable to store the number of GDCS genes were matched
                            count = 0
                            # As I want the count to be in the report before all the gene results, this string will
                            # store the specific sample information, and will be added to data once count is known
                            specific = str()
                            for gene in sorted(gdcs):
                                # As there are different genes present in the GDCS databases for each organism of
                                # interest, genes that did not match because they're absent in the specific database are
                                # indicated using an X
                                if gene not in [result for result in sample[analysistype].faidict]:
                                    specific += 'X,'
                                else:
                                    try:
                                        specific += '{p_id},'.format(p_id=sample[analysistype].blastresults[gene])
                                        # Report the necessary information for each gene result
                                        count += 1
                                    # If the gene was missing from the results attribute, add a - to the cell
                                    except (KeyError, AttributeError):
                                        specific += '-,'
                            # Determine whether the sample pass the necessary quality criteria:
                            # Pass, all GDCS, mean coverage greater than 20X coverage;
                            # ?: Indeterminate value;
                            # -: Fail value
                            # Allow one missing GDCS to still be considered a pass
                            if count >= len(sample[analysistype].faidict) - 1:
                                quality = '+'
                            else:
                                quality = '-'
                            # Add the count, mean depth with standard deviation, the pass/fail determination,
                            #  and the total number of GDCS genes as well as the results
                            data += '{hits}/{total},{fail},{gdcs}\n'\
                                .format(hits=str(count),
                                        total=len(sample[analysistype].faidict),
                                        fail=quality,
                                        gdcs=specific)
                        # Any samples with a best assembly of 'NA' are considered incomplete.
                        else:
                            data += '{},{},,,-\n'.format(sample.name, sample.general.closestrefseqgenus)
                    elif sample.general.closestrefseqgenus == 'NA':
                        data += '{}\n'.format(sample.name)
            # Write the header and data to file
            report.write(header)
            report.write(data)
        # Return the updated metadata object
        return metadata

    def sero_reporter(self, metadata, analysistype, reportpath):
        """
        Creates a report of the results
        """
        logging.info('Creating {} report'.format(analysistype))
        metadata = self.serotype_escherichia(metadata=metadata,
                                             analysistype=analysistype)
        # Create the path in which the reports are stored
        make_path(reportpath)
        header = 'Strain,Serotype\n'
        data = ''
        with open(os.path.join(reportpath, '{}.csv'.format(analysistype)), 'w') as report:
            for sample in metadata:
                if sample.general.bestassemblyfile != 'NA':
                    data += sample.name + ','
                    if sample[analysistype].blastresults:
                        # Set the O-type as either the appropriate attribute, or O-untypable
                        if ';'.join(sample.serosippr.o_set) == '-':
                            otype = 'O-untypeable'
                        else:
                            otype = '{oset} ({opid})'.format(oset=';'.join(sample.serosippr.o_set),
                                                             opid=sample.serosippr.best_o_pid)
                        # Same as above, but for the H-type
                        if ';'.join(sample.serosippr.h_set) == '-':
                            htype = 'H-untypeable'

                        else:
                            htype = '{hset} ({hpid})'.format(hset=';'.join(sample.serosippr.h_set),
                                                             hpid=sample.serosippr.best_h_pid)
                        serotype = '{otype}:{htype}'.format(otype=otype,
                                                            htype=htype)
                        # Populate the data string
                        data += serotype if serotype != 'O-untypeable:H-untypeable' else 'ND'
                        data += '\n'
                    else:
                        data += '\n'
            report.write(header)
            report.write(data)
        return metadata

    @staticmethod
    def serotype_escherichia(metadata, analysistype):
        """
        Create attributes storing the best results for the O and H types
        """
        for sample in metadata:
            # Initialise negative results to be overwritten when necessary
            sample[analysistype].best_o_pid = '-'
            sample[analysistype].o_genes = ['-']
            sample[analysistype].o_set = ['-']
            sample[analysistype].best_h_pid = '-'
            sample[analysistype].h_genes = ['-']
            sample[analysistype].h_set = ['-']
            if sample.general.bestassemblyfile != 'NA':
                if sample.general.closestrefseqgenus == 'Escherichia':
                    o = dict()
                    h = dict()
                    for result, percentid in sample[analysistype].blastresults.items():
                        if 'O' in result.split('_')[-1]:
                            o.update({result: float(percentid)})
                        if 'H' in result.split('_')[-1]:
                            h.update({result: float(percentid)})
                    # O
                    try:
                        sorted_o = sorted(o.items(), key=operator.itemgetter(1), reverse=True)
                        sample[analysistype].best_o_pid = str(sorted_o[0][1])

                        sample[analysistype].o_genes = [gene for gene, pid in o.items()
                                                        if str(pid) == sample[analysistype].best_o_pid]
                        sample[analysistype].o_set = \
                            list(set(gene.split('_')[-1] for gene in sample[analysistype].o_genes))
                    except (KeyError, IndexError):
                        pass
                    # H
                    try:
                        sorted_h = sorted(h.items(), key=operator.itemgetter(1), reverse=True)
                        sample[analysistype].best_h_pid = str(sorted_h[0][1])
                        sample[analysistype].h_genes = [gene for gene, pid in h.items()
                                                        if str(pid) == sample[analysistype].best_h_pid]
                        sample[analysistype].h_set = \
                            list(set(gene.split('_')[-1] for gene in sample[analysistype].h_genes))
                    except (KeyError, IndexError):
                        pass
        return metadata

    def alignprotein(self, sample, analysistype, target, targetfiles, records, program):
        """
        Create alignments of the sample nucleotide and amino acid sequences to the reference sequences
        :param sample: Metadata object
        :param analysistype: Current analysis type
        :param target: Current gene name
        :param targetfiles: List of all database files used in the analysis
        :param records: dictionary of Seq objects for all sequences in each database file
        :param program BLAST program used in the analyses
        :return: updated sample object
        """
        # Remove any gaps incorporated into the sequence
        sample[analysistype].targetsequence[target] = \
            sample[analysistype].targetsequence[target].replace('-', '')
        if program == 'blastn':
            # In order to properly translate the nucleotide sequence, BioPython requests that the sequence is a
            # multiple of three - not partial codons. Trim the sequence accordingly
            remainder = 0 - len(sample[analysistype].targetsequence[target]) % 3
            seq = sample[analysistype].targetsequence[target] if remainder == 0 \
                else sample[analysistype].targetsequence[target][:remainder]
            # Set the DNA and protein sequences of the target in the sample
            sample[analysistype].dnaseq[target] = Seq(seq, IUPAC.unambiguous_dna)
            # Translate the nucleotide sequence
            sample[analysistype].protseq[target] = str(sample[analysistype].dnaseq[target].translate())
        else:
            seq = sample[analysistype].targetsequence[target]
            sample[analysistype].protseq[target] = Seq(seq, IUPAC.protein)
        for targetfile in targetfiles:
            if program == 'blastn' or program == 'tblastn' or program == 'tblastx':
                # Trim the reference sequence to multiples of three
                refremainder = 0 - len(records[targetfile][target].seq) % 3
                refseq = str(records[targetfile][target].seq) if refremainder % 3 == 0 \
                    else str(records[targetfile][target].seq)[:refremainder]
                # Translate the nucleotide sequence of the reference sequence
                refdna = Seq(refseq, IUPAC.unambiguous_dna)
                try:
                    refprot = str(refdna.translate())
                # The ResFinder database seems to have at least one ambiguous DNA sequence, translate it appropriately
                except TranslationError:
                    # Translate the nucleotide sequence of the reference sequence
                    refdna = Seq(refseq, IUPAC.ambiguous_dna)
                    refprot = str(refdna.translate())
                # Use pairwise2 to perform a local alignment with the following parameters:
                # x     No match parameters. Identical characters have score of 1, otherwise 0.
                # s     Same open (-1)  and extend (-.1) gap penalties for both sequences
                ntalignments = pairwise2.align.localxs(seq, refseq, -1, -.1)
                # Use format_alignment to create a formatted alignment that is subsequently split on newlines e.g.
                '''
                ACCGT
                | ||
                A-CG-
                Score=3
                '''
                ntformat = (str(format_alignment(*ntalignments[0])).split('\n'))
                # Align the nucleotide sequence of the reference (ntalignments[2]) to the sample (ntalignments[0]).
                # If the corresponding bases match, add a |, otherwise a space
                ntalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                          zip(ntformat[0], ntformat[2])))
                # Create the nucleotide alignment: the sample sequence, the (mis)matches, and the reference sequence
                sample[analysistype].ntalign[target] = self.interleaveblastresults(ntformat[0], ntformat[2])
                # Regex to determine location of mismatches in the sequences
                count = 0
                sample[analysistype].ntindex[target] = str()
                for snp in re.finditer(' ', ntalignment):
                    # If there are many SNPs, then insert line breaks for every 10 SNPs
                    if count <= 10:
                        sample[analysistype].ntindex[target] += str(snp.start()) + ';'
                    else:
                        sample[analysistype].ntindex[target] += '\n' + str(snp.start()) + ';'
                        count = 0
                    count += 1
            else:
                refseq = str(records[targetfile][target].seq)
                # Translate the nucleotide sequence of the reference sequence
                refprot = Seq(refseq, IUPAC.protein)
            # Perform the same steps, except for the amino acid sequence
            aaalignments = pairwise2.align.localxs(sample[analysistype].protseq[target], refprot, -1, -.1)
            aaformat = (str(format_alignment(*aaalignments[0])).split('\n'))
            aaalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                      zip(aaformat[0], aaformat[2])))
            sample[analysistype].aaidentity[target] = '{:.2f}'\
                .format(float(aaalignment.count('|')) / float(len(aaalignment)) * 100)
            sample[analysistype].aaalign[target] = self.interleaveblastresults(aaformat[0], aaformat[2])
            count = 0
            sample[analysistype].aaindex[target] = str()
            for snp in re.finditer(' ', aaalignment):
                if count <= 10:
                    sample[analysistype].aaindex[target] += str(snp.start()) + ';'
                else:
                    sample[analysistype].aaindex[target] += '\n' + str(snp.start()) + ';'
                    count = 0
                count += 1
        return sample

    @staticmethod
    def interleaveblastresults(query, subject):
        """
        Creates an interleaved string that resembles BLAST sequence comparisons
        :param query: Query sequence
        :param subject: Subject sequence
        :return: Properly formatted BLAST-like sequence comparison
        """
        # Initialise strings to hold the matches, and the final BLAST-formatted string
        matchstring = str()
        blaststring = str()
        # Iterate through the query
        for i, bp in enumerate(query):
            # If the current base in the query is identical to the corresponding base in the reference, append a '|'
            # to the match string, otherwise, append a ' '
            if bp == subject[i]:
                matchstring += '|'
            else:
                matchstring += ' '
        # Set a variable to store the progress through the sequence
        prev = 0
        # Iterate through the query, from start to finish in steps of 60 bp
        for j in range(0, len(query), 60):
            # BLAST results string. The components are: current position (padded to four characters), 'OLC', query
            # sequence, \n, matches, \n, 'ref', subject sequence. Repeated until all the sequence data are present.
            """
            0000 OLC ATGAAGAAGATATTTGTAGCGGCTTTATTTGCTTTTGTTTCTGTTAATGCAATGGCAGCT
                     ||||||||||| ||| | |||| ||||||||| || ||||||||||||||||||||||||
                 ref ATGAAGAAGATGTTTATGGCGGTTTTATTTGCATTAGTTTCTGTTAATGCAATGGCAGCT
            0060 OLC GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
                     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                 ref GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
            """
            blaststring += '{} OLC {}\n         {}\n     ref {}\n' \
                .format('{:04d}'.format(j), query[prev:j + 60], matchstring[prev:j + 60], subject[prev:j + 60])
            # Update the progress variable
            prev = j + 60
        # Return the properly formatted string
        return blaststring

    @staticmethod
    def clean_object(metadata, analysistype):
        """
        Remove certain attributes from the object; they take up too much room on the .json report
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        """
        for sample in metadata:
            try:
                delattr(sample[analysistype], "targetnames")
            except AttributeError:
                pass
            try:
                delattr(sample[analysistype], "targets")
            except AttributeError:
                pass
            try:
                delattr(sample[analysistype], "dnaseq")
            except AttributeError:
                pass
            try:
                delattr(sample[analysistype], "protseq")
            except AttributeError:
                pass
