#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import run_subprocess, make_path, \
    combinetargets, MetadataObject, GenObject, printtime
from biotools.bbtools import kwargs_to_string
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastxCommandline, NcbiblastpCommandline, \
    NcbitblastnCommandline, NcbitblastxCommandline
from Bio.Application import ApplicationError
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from csv import DictReader
import multiprocessing
from glob import glob
import xlsxwriter
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

    def run_blast(self, metadata, analysistype, program, outfmt, evalue='1E-5', num_threads=12, num_alignments=1000000):
        """
        Runs BLAST on all the samples in the metadata object
        :param metadata: Metadata object
        :param analysistype: Name of analysis type
        :param program: BLAST program to use for the alignment
        :param outfmt; Custom fields to include in BLAST output
        :param evalue: e-value cut-off for BLAST analyses
        :param num_threads: Number of threads to use for BLAST analyses
        :param num_alignments: Number of alignments to perform in BLAST analyses
        :return: Updated metadata object
        """
        for sample in metadata:
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
                # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
                # remove the report, and run the blast analyses again
                if size == 0:
                    os.remove(sample[analysistype].report)
            except FileNotFoundError:
                pass
            # Split the extension from the file path
            db = os.path.splitext(sample[analysistype].combinedtargets)[0]
            # Create the command line argument using the appropriate BioPython BLAST wrapper
            if program == 'blastn':
                blast = self.blastn_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt)
            elif program == 'blastp':
                blast = self.blastp_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt)
            elif program == 'blastx':
                blast = self.blastx_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt)
            elif program == 'tblastn':
                blast = self.tblastn_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt)
            elif program == 'tblastx':
                blast = self.tblastx_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt)
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
                    raise
        # Return the updated metadata object
        return metadata

    @staticmethod
    def blastn_commandline(sample, analysistype, db, evalue, num_alignments, num_threads, outfmt):
        # BLAST command line call. Note the high number of default alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        blastn = NcbiblastnCommandline(query=sample.general.bestassemblyfile,
                                       db=db,
                                       evalue=evalue,
                                       num_alignments=num_alignments,
                                       num_threads=num_threads,
                                       outfmt=outfmt,
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
            # Open the sequence profile file as a dictionary
            blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
            resultdict = dict()
            # Initialise a dictionary to store all the target sequences
            sample[analysistype].targetsequence = dict()
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
            # Encountering the following error: # _csv.Error: field larger than field limit (131072)
            # According to https://stackoverflow.com/a/15063941, increasing the field limit should fix the issue
            csv.field_size_limit(sys.maxsize)
            # Open the sequence profile file as a dictionary
            blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
            # Initialise a dictionary to store all the target sequences
            sample[analysistype].targetsequence = dict()
            sample[analysistype].queryranges = dict()
            sample[analysistype].querypercent = dict()
            sample[analysistype].queryscore = dict()
            sample[analysistype].results = dict()
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
                target = row['subject_id']
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
                        # seen before, append it. Otherwise, update the previous ranges with the new, longer range as
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
                            # Update the lsit with the blast results
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
        extended = False
        for sample in metadata:
            # Reset the column to zero
            col = 0
            # Initialise a list to store all the data for each strain
            data = list()
            # Initialise a list of all the headers with 'Strain'
            headers = ['Strain']
            if sample[analysistype].targetnames != 'NA':
                # Append the sample name to the data list only if the script could find targets
                data.append(sample.name)
                if sample[analysistype].blastresults != 'NA':
                    for target in sorted(sample[analysistype].targetnames):
                        # Add the name of the gene to the header
                        headers.append(target)
                        try:
                            # Append the percent identity to the data list
                            data.append(str(sample[analysistype].blastresults[target]))
                            # Only if the alignment option is selected, for inexact results, add alignments
                            if align and sample[analysistype].blastresults[target] != 100.00:
                                # Align the protein (and nucleotide) sequences to the reference
                                sample = self.alignprotein(sample, analysistype, target, targetfiles, records, program)
                                if not extended:
                                    if program == 'blastn':
                                        # Add the appropriate headers
                                        headers.extend(['aa_Sequence',
                                                        'aa_Alignment',
                                                        'aa_SNP_location',
                                                        'nt_Alignment',
                                                        'nt_SNP_location'
                                                        ])
                                    else:
                                        headers.extend(['aa_Sequence',
                                                        'aa_Alignment',
                                                        'aa_SNP_location',
                                                        ])
                                    extended = True
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
                                    data.extend([record.format('fasta'),
                                                 sample[analysistype].aaalign[target],
                                                 sample[analysistype].aaindex[target],
                                                 sample[analysistype].ntalign[target],
                                                 sample[analysistype].ntindex[target]
                                                 ])
                                else:
                                    data.extend([record.format('fasta'),
                                                 sample[analysistype].aaalign[target],
                                                 sample[analysistype].aaindex[target],
                                                 ])
                        # If there are no blast results for the target, add a '-'
                        except (KeyError, TypeError):
                            data.append('-')
                # If there are no blast results at all, add a '-'
                else:
                    data.append('-')
            # Write the header to the spreadsheet
            for header in headers:
                worksheet.write(row, col, header, bold)
                # Set the column width based on the longest header
                try:
                    columnwidth[col] = len(header)if len(header) > columnwidth[col] else columnwidth[col]
                except KeyError:
                    columnwidth[col] = len(header)
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
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
                worksheet.set_row(row, max(totallines) * 12)
            else:
                worksheet.set_row(row, 1)
            # Increase the row counter for the next strain's data
            row += 1
        # Close the workbook
        workbook.close()
        # Return the updated metadata object
        return metadata

    def resfinder_reporter(self, metadata, analysistype, targetfolders, reportpath, align, targetfiles, records,
                           program):
        """
        Custom reports for ResFinder analyses. These reports link the gene(s) found to their resistance phenotypes
        :param metadata: Metadata object
        :param analysistype: Current analysis type
        :param targetfolders: List of all folders with targets used in the analyses
        :param reportpath: Path of folder in which report is to be created
        :param align: Boolean of whether alignments between query and subject sequences are desired
        :param targetfiles: List of all files used in the analyses
        :param records: Dictionary of SeqIO parsed sequence records
        :param program: BLAST program used in the analyses
        :return: Updated metadata object
        """
        target_dir = str()
        for folder in targetfolders:
            target_dir = folder
        genedict, altgenedict = ResistanceNotes.notes(target_dir)
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
            sample[analysistype].pipelineresults = list()
            # Process the sample only if the script could find targets
            if sample[analysistype].blastlist != 'NA' and sample[analysistype].blastlist:
                for result in sample[analysistype].blastlist:
                    # Set the name to avoid writing out the dictionary[key] multiple times
                    name = result['subject_id']
                    # Use the ResistanceNotes gene name extraction method to get the necessary variables
                    gname, genename, accession, allele = ResistanceNotes.gene_name(name)
                    # Initialise a list to store all the data for each strain
                    data = list()
                    # Determine the name of the gene to use in the report and the resistance using the resistance
                    # method
                    finalgene, resistance = ResistanceNotes.resistance(gname, genename, genedict, altgenedict)
                    # Append the necessary values to the data list
                    data.append(finalgene)
                    data.append(allele)
                    data.append(resistance)
                    percentid = result['percentidentity']
                    data.append(percentid)
                    data.append(result['alignment_fraction'])
                    data.append(result['query_id'])
                    data.append('...'.join([str(result['low']), str(result['high'])]))
                    # Populate the .pipelineresults attribute for compatibility with the assembly pipeline
                    sample[analysistype].pipelineresults.append(
                        '{rgene} ({pid}%) {rclass}'.format(rgene=finalgene,
                                                           pid=percentid,
                                                           rclass=resistance))
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
            sample[analysistype].protseq[target] = Seq(seq,
                                                       IUPAC.protein)
        for targetfile in targetfiles:
            if program == 'blastn' or program == 'tblastn' or program == 'tblastx':
                # Trim the reference sequence to multiples of three
                refremainder = 0 - len(records[targetfile][target].seq) % 3
                refseq = str(records[targetfile][target].seq) if refremainder % 3 == 0 \
                    else str(records[targetfile][target].seq)[:refremainder]
                # Translate the nucleotide sequence of the reference sequence
                refdna = Seq(refseq, IUPAC.unambiguous_dna)
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
            except KeyError:
                pass
            try:
                delattr(sample[analysistype], "targets")
            except KeyError:
                pass
            try:
                delattr(sample[analysistype], "dnaseq")
            except KeyError:
                pass
            try:
                delattr(sample[analysistype], "protseq")
            except KeyError:
                pass


class Parser(object):

    def main(self):
        """
        Run the parsing methods
        """
        self.target_find()
        self.strainer()
        self.metadata_populate()

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
        assert self.targets, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                             'to ensure that your target path is correct'.format(self.targetpath)

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

    def metadata_populate(self):
        """
        Populate the :analysistype GenObject
        """
        for metadata in self.metadata:
            # Create and populate the :analysistype attribute
            setattr(metadata, self.analysistype, GenObject())

            metadata[self.analysistype].targets = self.targets
            metadata[self.analysistype].combinedtargets = self.combinedtargets
            metadata[self.analysistype].targetpath = self.targetpath
            metadata[self.analysistype].targetnames = sequencenames(self.combinedtargets)
            try:
                metadata[self.analysistype].reportdir = os.path.join(metadata.general.outputdirectory,
                                                                     self.analysistype)
            except (AttributeError, KeyError):
                metadata[self.analysistype].reportdir = self.reportpath

    def __init__(self, args):
        self.analysistype = args.analysistype
        self.sequencepath = os.path.join(args.sequencepath)
        self.targetpath = os.path.join(args.targetpath)
        if not os.path.isdir(self.targetpath):
            self.targetpath = self.targetpath.split('_')[0]
        assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}' \
            .format(self.targetpath)
        self.reportpath = os.path.join(args.reportpath)
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


class ResistanceNotes(object):

    @staticmethod
    def notes(targetpath):
        """
        Populates resistance dictionary with resistance class: base gene name
        :param targetpath: Directory in which the notes.txt file is located
        :return: the resistance dictionaries
        """
        # Create a set of all the resistance classes and the base names
        genedict = dict()
        altgenedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            res_class = str()
            for line in notes:
                # Create entries for each class - these are on lines beginning with '#' e.g. #Rifampicin resistance
                if line.startswith('#'):
                    res_class = line.split(' resistance')[0].lstrip('#').replace(':', '').rstrip()
                    # Initialise the dictionary as a set for the resistance class
                    genedict[res_class] = set()
                else:
                    # aac(6')-III:Aminoglycoside resistance: yields a base gene name of aac
                    gene = line.split(':')[0].split('(')[0].split('-')[0].rstrip()
                    # Aminoglycoside resistance
                    if 'Alternate name' in line:
                        # Extract the full gene name e.g. aac(6')-III:Aminoglycoside resistance: yields aac(6')-III
                        full_gene = line.split(':')[0]
                        # There are two formats for the alternate name in the file. Account for both
                        try:
                            alternate = line.split(';')[1].rstrip().lstrip()
                        except IndexError:
                            alternate = line.split('name ')[1].rstrip().lstrip()
                        # Populate the dictionaries
                        genedict[res_class].add(alternate)
                        altgenedict[full_gene] = alternate
                    else:
                        genedict[res_class].add(gene)
        return genedict, altgenedict

    @staticmethod
    def gene_name(name):
        """
        Split the FASTA header string into its components, including gene name, allele, and accession
        :param name: FASTA header
        :return: gname, genename, accession, allele: name of gene. Often the same as genename, but for certain entries
        it is longer, full gene name, accession, and allele extracted from the FASTA header
        """
        if 'Van' in name or 'mcr' in name or 'aph' in name or 'ddlA' in name or 'ant' in name:
            try:
                if name == "ant(3'')_Ih_aac(6')_IId_1_AF453998":
                    # >aac(3)_Ib_aac(6')_Ib_1_AF355189 yields gname, genename: aac(3)-Ib-aac(6')-Ib, allele:1,
                    # accession: AF355189
                    gene1, version1, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}-{g2}-{v2}'.format(g1=gene1,
                                                         v1=version1,
                                                         g2=gene2,
                                                         v2=version2)
                    genename = gname
                elif name == 'ant(3'')_Ia_1_X02340':
                    # >ant(3'')_Ia_1_X02340
                    gene, version, allele, accession = name.split('_')
                    gname = '{g}-{v}'.format(g=gene,
                                             v=version)
                    genename = gname
                elif 'mcr_3' in name or 'mcr_2' in name or 'mcr_1.10' in name:
                    # >mcr_3.3_1_NG055492 yields genename, gname: mcr-3, allele: 1, accession: NG055492
                    gene, combinedversion, allele, accession = name.split('_')
                    version = combinedversion.split('.')[0]
                    gname = '{gene}-{version}'.format(gene=gene,
                                                      version=version)
                    genename = gname

                else:
                    # Allow for an additional part to the gene name aph(3'')_Ib_5_AF321551 yields gname: aph(3''),
                    # genename: aph(3'')-Ib, allele: 5, accession AF321551
                    try:
                        pregene, postgene, allele, accession = name.split('_')
                        gname = '{pre}-{post}'.format(pre=pregene,
                                                      post=postgene)
                        genename = gname
                    except ValueError:
                        # Allow for underscores in the accession: aac(2')_Ie_1_NC_011896 yields gname: aac(2'),
                        # genename:  aac('2)-1e, allele: 1, accession NC_011896
                        pregene, postgene, allele, preaccession, postaccession = name.split('_')
                        genename = '{pre}-{post}'.format(pre=pregene,
                                                         post=postgene)
                        accession = '{pre}_{post}'.format(pre=preaccession,
                                                          post=postaccession)
                        gname = pregene
            except ValueError:
                # VanC_2_DQ022190
                genename, allele, accession = name.split('_')
                gname = genename
        else:
            if 'bla' in name or 'aac' in name or 'ARR' in name or 'POM' in name:
                if 'OKP' in name or 'CTX' in name or 'OXY' in name:
                    # >blaOKP_B_11_1_AM051161 yields gname: blaOKP-B-11, genename: blaOXP, allele: 1,
                    # accession: AM051161
                    gene, version1, version2, allele, accession = name.split('_')
                    gname = '{g}-{v1}-{v2}'.format(g=gene,
                                                   v1=version1,
                                                   v2=version2)
                    genename = gname
                elif 'CMY' in name:
                    # >blaCMY_12_1_Y16785 yields gname, genename: blaCMY, allele: 12
                    gname, allele, version, accession = name.split('_')
                    genename = gname
                elif name == "aac(3)_Ib_aac(6')_Ib_1_AF355189":
                    # >aac(3)_Ib_aac(6')_Ib_1_AF355189 yields gname, genename: aac(3)-Ib-aac(6')-Ib, allele:1,
                    # accession: AF355189
                    gene1, version1, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}-{g2}-{v2}'.format(g1=gene1,
                                                         v1=version1,
                                                         g2=gene2,
                                                         v2=version2)
                    genename = gname
                elif 'alias' in name:
                    # >blaSHV_5a_alias_blaSHV_9_1_S82452
                    gene1, version1, alias, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}'.format(g1=gene1,
                                               v1=version1)
                    genename = gname
                else:
                    # Split the name on '_'s: ARR-2_1_HQ141279; gname, genename: ARR-2, allele: 1, accession: HQ141279
                    try:
                        genename, allele, accession = name.split('_')
                        gname = genename
                    except ValueError:
                        try:
                            # >blaACC_1_2_AM939420 yields gname: blaACC-1, genename: blaACC, allele: 2,
                            # accession: AM939420
                            genename, version, allele, accession = name.split('_')
                            gname = '{g}-{v}'.format(g=genename,
                                                     v=version)
                        except ValueError:
                            # >aac(2')_Ie_1_NC_011896 yields gname, genename: aac(2')-Ie, allele: 1,
                            # accession: NC_011896
                            genename, version, allele, preaccession, postaccession = name.split('_')
                            gname = '{g}-{v}'.format(g=genename,
                                                     v=version)
                            genename = gname
                            accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                          postaccess=postaccession)
            else:
                # Split the name on '_'s: ARR-2_1_HQ141279; gname, genename: ARR-2, allele: 1, accession: HQ141279
                try:
                    genename, allele, accession = name.split('_')
                    gname = genename
                # Some names have a slightly different naming scheme:
                except ValueError:
                    # tet(44)_1_NZ_ABDU01000081 yields gname, genename: tet(44), allele: 1,
                    # accession: NZ_ABDU01000081
                    genename, allele, preaccession, postaccession = name.split('_')
                    accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                  postaccess=postaccession)
                    gname = genename
        return gname, genename, accession, allele

    @staticmethod
    def resistance(gname, genename, genedict, altgenedict):
        """
        Extracts the resistance phenotype from the dictionaries using the gene name
        :param gname: Name of gene. Often the same as genename, but for certain entries it is longer
        e.g. blaOKP-B-15 instead of blaOKP
        :param genename: Name of gene e.g. blaOKP
        :param genedict: Dictionary of gene:resistance
        :param altgenedict: Dictionary of gene alternate name:resistance
        :return: finalgene, finalresistance: gene name and associated resistance phenotype to be used in the report
        """
        # Use the same string splitting method as above to re-create the base name of the gene from the gene name
        # e.g. aph(3'')-Ib yields aph
        gene = gname.split('(')[0].split('-')[0]
        # Set the final gene name as the genename variable - unless there is an alternative name (see below)
        finalgene = genename
        # Initialise the finalresistance variable
        finalresistance = str()
        # Iterate through all the resistance classes, and gene sets associated with each resistance class
        for resistance, gene_set in genedict.items():
            # Determine if the gene is present in the resistance class-associated gene set
            if gene in gene_set or genename in gene_set:
                # Set the resistance to return as the current resistance class
                finalresistance = resistance
                # Determine if there is an alternative name for this gene in the notes.txt file
                try:
                    alt_gene = altgenedict[gname]
                    # Set the final name of the gene to be gene name (alternative name)
                    finalgene = '{namegene} ({genealt})'.format(namegene=genename,
                                                                genealt=alt_gene)
                # Otherwise use the gene name
                except KeyError:
                    finalgene = genename
        return finalgene, finalresistance


def sequencenames(contigsfile):
    """
    Takes a multifasta file and returns a list of sequence names
    :param contigsfile: multifasta of all sequences
    :return: list of all sequence names
    """
    sequences = list()
    for record in SeqIO.parse(open(contigsfile, "rU", encoding="iso-8859-15"), "fasta"):
        sequences.append(record.id)
    return sequences


def objector(kw_dict, start):
    metadata = MetadataObject()
    for key, value in kw_dict.items():
        setattr(metadata, key, value)
    # Set the analysis type based on the arguments provided
    if metadata.resfinder is True:
        metadata.analysistype = 'resfinder'
    elif metadata.virulencefinder is True:
        metadata.analysistype = 'virulence'
    # Warn that only one type of analysis can be perfomed at a time
    elif metadata.resfinder is True and metadata.virulencefinder is True:
        printtime('Cannot perform ResFinder and VirulenceFinder simultaneously. Please choose only one '
                  'of the -R and -v flags', start)
    # Default to GeneSeekr
    else:
        metadata.analysistype = 'geneseekr'
    # Add the start time variable to the object
    metadata.start = start
    return metadata


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
        printtime('{at} analyses complete'.format(at=self.analysistype), self.start)

    def blast_db(self):
        """
        Make blast databases (if necessary)
        """
        printtime('Creating {at} blast databases as required'
                  .format(at=self.analysistype),
                  self.start)
        self.geneseekr.makeblastdb(self.combinedtargets,
                                   self.program)
        self.targetfolders, self.targetfiles, self.records = \
            self.geneseekr.target_folders(self.metadata,
                                          self.analysistype)

    def run_blast(self):
        """
        Perform BLAST analyses
        """
        printtime('Performing {program} analyses on {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
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
        printtime('Parsing {program} results for {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
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
        printtime('Creating {at} reports'.format(at=self.analysistype), self.start)
        if 'resfinder' in self.analysistype:
            # ResFinder-specific report
            self.metadata = self.geneseekr.resfinder_reporter(self.metadata,
                                                              self.analysistype,
                                                              self.targetfolders,
                                                              self.reportpath,
                                                              self.align,
                                                              self.targetfiles,
                                                              self.records,
                                                              self.program)
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

    def __init__(self, args, analysistype='geneseekr', cutoff=70, program='blastn'):
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
        except (AttributeError, KeyError):
            self.targetpath = args.targetpath
        self.reportpath = args.reportpath
        try:
            self.metadata = args.runmetadata.samples
            parse = Parser(self)
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
        printtime('Performing {program} analyses on {at} targets'
                  .format(program=self.program,
                          at=self.analysistype),
                  self.start)
