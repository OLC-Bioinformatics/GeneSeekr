#!/usr/bin/env python
import os
import sys
import time
import signal
from collections import defaultdict
from Bio.Application import _Option, AbstractCommandline, _Switch
from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing import Pool
__author__ = 'mike knowles'

__doc__ = 'The purpose of this set of modules is to improve upon earlier development of ARMISeekr.py and eventually' \
          'to include generalized functionality for with OOP for GeneSeekr'

class KeyboardInterruptError(Exception): pass


class MakeBlastDB(AbstractCommandline):
    """Base makeblastdb wrapper"""
    def __init__(self, cmd='makeblastdb', **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file prefix for db.",
                    filename=True,
                    equate=False),
            _Option(["-in", "db"],
                    "The sequence create db with.",
                    filename=True,
                    equate=False),  # Should this be required?
            _Option(["-dbtype", "dbtype"],
                    "Molecule type of target db (string, 'nucl' or 'prot').",
                    equate=False)]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):  # deal with mangled names
        cls_name = cls.__name__.lstrip('_')
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class ARMISeekr(object):

    def yeah(self, reset=None):
        """
        :type reset: int
        :rtype: yeah
        """
        if reset is not None:
            self.count = 1
        if self.count == 1:
            sys.stdout.write('[{}] 1 ( \xE2\x80\xA2_\xE2\x80\xA2)'.format(time.strftime("%H:%M:%S")))
        elif self. count % 3 == 0:
            sys.stdout.write('\r[{}] {} (\xE2\x8C\x90\xE2\x96\xA0_\xE2\x96\xA0) #Yeeeaaahhhh'
                             .format(time.strftime("%H:%M:%S"), self.count))
        elif self.count % 2 == 0:
            sys.stdout.write('\r[{}] {} ( \xE2\x80\xA2_\xE2\x80\xA2)>\xE2\x8C\x90\xE2\x96\xA0-\xE2\x96\xA0'
                             .format(time.strftime("%H:%M:%S"), self.count))
        else:
            sys.stdout.write('\r[{}] {} ( \xE2\x80\xA2_\xE2\x80\xA2)'.format(time.strftime("%H:%M:%S"), self.count))
        self.count += 1

    def __init__(self, subject, query, threads=12):
        """:type subject: list of genes
           :type query: list of target genomes"""
        assert isinstance(subject, list), 'Subject is not a list "{0!r:s}"'.format(subject)
        assert isinstance(query, list), 'Query is not a list"{0!r:s}"'.format(query)
        self.count, self.subject, self.query, self.threads = 0, subject, query, threads
        self.cutoff, self.genelist = 70, []
        self.db = map((lambda x: os.path.splitext(x)[0]), subject)  # remove the file extension for easier globing
        self.plus = dict((target, defaultdict(list)) for target in self.query)  # Initialize :return dict
        print '[{}] GeneSeekr input is path with {} files'.format(time.strftime("%H:%M:%S"), len(query))
        print "[{}] Creating necessary databases for BLAST".format(time.strftime("%H:%M:%S"))
        pool = Pool(self.threads)
        try:
            pool.map(makeblastdb, zip(self.subject, self.db))
        except KeyboardInterrupt:
            print "[{0:s}] Got ^C while pool mapping, terminating the pool".format(time.strftime("%H:%M:%S"))
            pool.terminate()
            print 'pool is terminated'
            sys.exit(127)
        except Exception, e:
            print "[{0:s}] Got exception: {1!r:s}, terminating the pool".format(e, time.strftime("%H:%M:%S"))
            pool.terminate()
            print "[{0:s}] Pool is terminated".format(time.strftime("%H:%M:%S"))
            sys.exit(127)
        print "\r[{0}] BLAST database(s) created".format(time.strftime("%H:%M:%S"))

    def _blast(self, (fasta, db)):
        try:
            blastn = NcbiblastnCommandline(query=fasta,
                                           db=db,
                                           evalue=10,
                                           outfmt="'6 sseqid nident slen'",
                                           perc_identity=self.cutoff)
            stdout, stderr = blastn()
            if stdout != '':
                return [[fasta, aln[0][4:], abs(float(aln[1]) / float(aln[2]))]
                        for aln in [hsp.split('\t')
                                    for hsp in stdout.rstrip().split("\n")]
                        if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]
        except KeyboardInterrupt:
            raise KeyboardInterruptError()

    def mpblast(self, cutoff=70):
        assert isinstance(cutoff, int), u'Cutoff is not an integer {0!r:s}'.format(cutoff)
        self.cutoff = cutoff
        print "[{0:s}] Now performing and parsing BLAST database searches".format(time.strftime("%H:%M:%S"))
        start = time.time()
        p = Pool(self.threads)
        for genes in self.db:
            try:
                mapblast = p.map(self._blast, [(genome, genes) for genome in self.query])
                for fastaline in mapblast:
                    if fastaline is not None:  # if the returned list contains [genome, gene, value]
                        for fasta, gene, v in fastaline:  # unpack
                            if gene not in self.genelist:
                                self.genelist.append(gene)  # create list of all genes in anaylsis
                            self.plus[fasta][gene].append(v)
                            self.plus[fasta][gene].sort()
            except KeyboardInterrupt:
                print "[{0:s}] Got ^C while pool mapping, terminating the pool".format(time.strftime("%H:%M:%S"))
                p.terminate()
                print 'pool is terminated'
                sys.exit(127)
            except Exception, e:
                print "[{0:s}] Got exception: {1!r:s}, terminating the pool".format(e, time.strftime("%H:%M:%S"))
                p.terminate()
                print "[{0:s}] Pool is terminated".format(time.strftime("%H:%M:%S"))
                sys.exit(127)

        print "[{}] Now compiling BLAST database results".format(time.strftime("%H:%M:%S"))
        end = time.time() - start
        print "[{0:s}] Elapsed time for GeneSeekr is {1:0d}m {2:0d}s with {3:0.2f}s per genome".format(
            time.strftime("%H:%M:%S"), int(end) / 60, int(end) % 60, end / float(len(self.query)))
        return self.plus

    def csvwriter(self, out, name):
        assert isinstance(out, str), u'Output location is not a string {0!r:s}'.format(out)
        assert isinstance(name, str), u'Output name is not a string {0!r:s}'.format(name)
        assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
        print "[{}] Writing CSV and JSON to output directory".format(time.strftime("%H:%M:%S"))
        self.genelist.sort()
        rowcount, row = 0, 'Strain,'
        row += ', '.join(self.genelist)
        for genomerow in sorted(self.plus):
            row += '\n{}'.format(os.path.split(os.path.splitext(genomerow)[0])[1].replace('_filteredAssembled', ""))
            for genename in self.genelist:
                row += ',' + (lambda x, y: ' '.join(map(str, x[y])) if y in x else 'N')(self.plus[genomerow], genename)
                # Add the allele numbers to the row for the appropriate gene, otherwise return N
        with open("%s/%s_results_%s.csv" % (out, name, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
            csvfile.write(row)


def makeblastdb((fasta, db)):
    try:
        if not os.path.isfile('{}.nhr'.format(db)):  # add nhr for searching
            assert os.path.isfile(fasta)  # check that the fasta has been specified properly
            MakeBlastDB(db=fasta, out=db, dbtype='nucl')()  # Use MakeBlastDB above
        return 0
    except KeyboardInterrupt:
            raise KeyboardInterruptError()
