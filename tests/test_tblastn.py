#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.geneseekr import GeneSeekr
import geneseekr.tblastn as tblastn
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
    v.targetpath = os.path.join(datapath, 'databases', 'resfinder')
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


tblastn_method = method_init(variables(), 'resfinder', 'tblastn', True, True, tblastn.tBLASTn)


def test_parser():
    assert os.path.basename(tblastn_method.targets[0]) == 'beta-lactam.tfa'


def test_combined_files():
    assert os.path.isfile(tblastn_method.combinedtargets)


def test_strains():
    assert os.path.isfile(tblastn_method.strains[0])


def test_strain():
    assert os.path.basename(tblastn_method.strains[0]) == 'amr_test.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(tblastn_method.combinedtargets)
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.nsq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(tblastn_method.metadata,
                                 tblastn_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'resfinder'


def test_targetfiles():
    assert targetfiles[0] == tblastn_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['ampH_2_HQ586946']


def test_tblastn():
    tblastn_method.metadata = geneseekr.run_blast(tblastn_method.metadata,
                                                  tblastn_method.analysistype,
                                                  tblastn_method.program,
                                                  tblastn_method.outfmt,
                                                  evalue=tblastn_method.evalue,
                                                  num_threads=tblastn_method.cpus)


def test_tblastn_report(variables):
    global tblastn_report
    tblastn_report = os.path.join(variables.reportpath, 'amr_test_tblastn.csv')
    assert os.path.isfile(tblastn_report)


def test_tblastn_results():
    with open(tblastn_report) as blast_results:
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    tblastn_method.metadata = geneseekr.unique_parse_blast(tblastn_method.metadata,
                                                           tblastn_method.analysistype,
                                                           tblastn_method.fieldnames,
                                                           tblastn_method.cutoff,
                                                           tblastn_method.program)
    for sample in tblastn_method.metadata:
        assert sample.resfinder.queryranges['contig2'] == [[1, 264]]


def test_filter():
    tblastn_method.metadata = geneseekr.filter_unique(tblastn_method.metadata,
                                                      tblastn_method.analysistype)
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    tblastn_method.metadata = geneseekr.dict_initialise(tblastn_method.metadata,
                                                        tblastn_method.analysistype)
    for sample in tblastn_method.metadata:
        assert type(sample.resfinder.protseq) is dict


def test_report_creation():
    tblastn_method.metadata = geneseekr.resfinder_reporter(tblastn_method.metadata,
                                                           tblastn_method.analysistype,
                                                           targetfolders,
                                                           tblastn_method.reportpath,
                                                           tblastn_method.align,
                                                           tblastn_method.targetfiles,
                                                           tblastn_method.records,
                                                           tblastn_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(tblastn_method.reportpath, 'resfinder_tblastn.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_report_row():
    for sample in tblastn_method.metadata:
        assert sorted(sample.resfinder.sampledata)[0] == \
               ['blaOXA', '1', 'Beta-Lactamase', 94.34, 99.62, 'contig2', '1...264', '-']


def test_parse_results():
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastresults['blaOXA_427_1_KX827604'] == 94.34


def test_aaseq():
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_combined_targets_clean():
    os.remove(tblastn_method.combinedtargets)


def test_makeblastdb_clean(variables):
    databasefiles = glob(os.path.join(variables.targetpath, 'combinedtargets.n*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_tblastn_report():
    os.remove(tblastn_report)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)
