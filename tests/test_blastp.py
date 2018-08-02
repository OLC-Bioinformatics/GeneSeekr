#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.geneseekr import GeneSeekr
import geneseekr.blastp as blastp
import multiprocessing
from glob import glob
from time import time
import pytest
import os

test_path = os.path.abspath(os.path.dirname(__file__))

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = MetadataObject()
    datapath = os.path.join(test_path, 'testdata')
    v.sequencepath = os.path.join(datapath, 'aa_sequences')
    v.targetpath = os.path.join(datapath, 'databases', 'card_aa')
    v.reportpath = os.path.join(datapath, 'reports')
    v.cutoff = 70
    v.evalue = '1E-05'
    v.align = False
    v.unique = False
    v.resfinder = False
    v.virulencefinder = False
    v.numthreads = multiprocessing.cpu_count()
    v.start = time()
    return v


def variable_update():
    global method
    method = method_init(variables())


@pytest.fixture()
def method_init(variables, analysistype, program, align, unique, dest):
    global method
    variables.analysistype = analysistype
    variables.program = program
    variables.align = align
    variables.unique = unique
    method = dest(variables)
    return method


blastp_method = method_init(variables(), 'geneseekr', 'blastp', True, True, blastp.BLASTp)


def test_parser():
    assert os.path.basename(blastp_method.targets[0]) == 'amr.tfa'


def test_combined_files():
    assert os.path.isfile(blastp_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastp_method.strains[0])


def test_strain():
    assert os.path.basename(blastp_method.strains[0]) == 'amr_test.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(blastp_method.combinedtargets,
                          dbtype='prot')
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.psq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(blastp_method.metadata,
                                 blastp_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'card_aa'


def test_targetfiles():
    assert targetfiles[0] == blastp_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['yojI']


def test_blastp():
    blastp_method.metadata = geneseekr.run_blast(blastp_method.metadata,
                                                 blastp_method.analysistype,
                                                 blastp_method.program,
                                                 blastp_method.outfmt,
                                                 evalue=blastp_method.evalue,
                                                 num_threads=blastp_method.cpus)


def test_blastp_report(variables):
    global blastp_report
    blastp_report = os.path.join(variables.reportpath, 'amr_test_blastp.csv')
    assert os.path.isfile(blastp_report)


def test_blastp_results():
    with open(blastp_report) as blast_results:
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    blastp_method.metadata = geneseekr.unique_parse_blast(blastp_method.metadata,
                                                          blastp_method.analysistype,
                                                          blastp_method.fieldnames,
                                                          blastp_method.cutoff,
                                                          blastp_method.program)
    for sample in blastp_method.metadata:
        assert sample.geneseekr.queryranges['contig1'] == [[1, 547]]


def test_filter():
    blastp_method.metadata = geneseekr.filter_unique(blastp_method.metadata,
                                                     blastp_method.analysistype)
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    blastp_method.metadata = geneseekr.dict_initialise(blastp_method.metadata,
                                                       blastp_method.analysistype)
    for sample in blastp_method.metadata:
        assert type(sample.geneseekr.protseq) is dict


def test_report_creation():
    blastp_method.metadata = geneseekr.reporter(blastp_method.metadata,
                                                blastp_method.analysistype,
                                                blastp_method.reportpath,
                                                blastp_method.align,
                                                blastp_method.targetfiles,
                                                blastp_method.records,
                                                blastp_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(blastp_method.reportpath, 'geneseekr_blastp.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_parse_results():
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastresults['OXA_12'] == 94.19


def test_aaseq():
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastlist[0]['query_sequence'][:4] == 'MELL' or \
               sample.geneseekr.blastlist[0]['query_sequence'][:4] == 'MSRI'


def test_combined_targets_clean():
    os.remove(blastp_method.combinedtargets)


def test_makeblastdb_clean(variables):
    databasefiles = glob(os.path.join(variables.targetpath, 'combinedtargets.p*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_blastp_report():
    os.remove(blastp_report)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)


def test_remove_report_path():
    os.rmdir(blastp_method.reportpath)
