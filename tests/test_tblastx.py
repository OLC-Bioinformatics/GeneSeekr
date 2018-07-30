#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
import multiprocessing
from glob import glob
from time import time
import pytest
import sys
import os

test_path = os.path.abspath(os.path.dirname(__file__))
script_path = os.path.join(test_path, '..')
sys.path.append(script_path)
method_path = os.path.join(script_path, 'methods')
sys.path.append(method_path)
bin_path = os.path.join(script_path, 'bin')
sys.path.append(bin_path)

from methods.geneseekr import GeneSeekr, Parser, sequencenames
import tblastx

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = MetadataObject()
    datapath = os.path.join(test_path, 'testdata')
    v.sequencepath = os.path.join(datapath, 'sequences')
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


tblastx_method = method_init(variables(), 'resfinder', 'tblastx', True, True, tblastx.tBLASTx)


def test_parser():
    assert os.path.basename(tblastx_method.targets[0]) == 'beta-lactam.tfa'


def test_combined_files():
    assert os.path.isfile(tblastx_method.combinedtargets)


def test_strains():
    assert os.path.isfile(tblastx_method.strains[0])


def test_strain():
    assert os.path.basename(tblastx_method.strains[0]) == '2018-SEQ-0552.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(tblastx_method.combinedtargets)
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.nsq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(tblastx_method.metadata,
                                 tblastx_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'resfinder'


def test_targetfiles():
    assert targetfiles[0] == tblastx_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['ampH_2_HQ586946']


def test_tblastx():
    tblastx_method.metadata = geneseekr.run_blast(tblastx_method.metadata,
                                                  tblastx_method.analysistype,
                                                  tblastx_method.program,
                                                  tblastx_method.outfmt,
                                                  evalue=tblastx_method.evalue,
                                                  num_threads=tblastx_method.cpus)


def test_tblastx_report(variables):
    global tblastx_report
    tblastx_report = os.path.join(variables.reportpath, '2018-SEQ-0552_tblastx.csv')
    assert os.path.isfile(tblastx_report)


def test_tblastx_results():
    with open(tblastx_report) as blast_results:
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    tblastx_method.metadata = geneseekr.unique_parse_blast(tblastx_method.metadata,
                                                           tblastx_method.analysistype,
                                                           tblastx_method.fieldnames,
                                                           tblastx_method.cutoff,
                                                           tblastx_method.program)
    for sample in tblastx_method.metadata:
        assert sample.resfinder.queryranges['Contig_54_76.3617'] == [[11054, 11848]]


def test_filter():
    tblastx_method.metadata = geneseekr.filter_unique(tblastx_method.metadata,
                                                      tblastx_method.analysistype)
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    tblastx_method.metadata = geneseekr.dict_initialise(tblastx_method.metadata,
                                                        tblastx_method.analysistype)
    for sample in tblastx_method.metadata:
        assert type(sample.resfinder.protseq) is dict


def test_report_creation():
    tblastx_method.metadata = geneseekr.resfinder_reporter(tblastx_method.metadata,
                                                           tblastx_method.analysistype,
                                                           targetfolders,
                                                           tblastx_method.reportpath,
                                                           tblastx_method.align,
                                                           tblastx_method.targetfiles,
                                                           tblastx_method.records,
                                                           tblastx_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(tblastx_method.reportpath, 'resfinder_tblastx.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_report_row():
    for sample in tblastx_method.metadata:
        assert sorted(sample.resfinder.sampledata)[0] == \
               ['ampH', '2', 'Beta-Lactamase', 93.96, 100.0, 'Contig_54_76.3617', '11054...11848', '-']


def test_parse_results():
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastresults['ampH_2_HQ586946'] == 93.96


def test_aaseq():
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_combined_targets_clean():
    os.remove(tblastx_method.combinedtargets)


def test_makeblastdb_clean(variables):
    databasefiles = glob(os.path.join(variables.targetpath, 'combinedtargets.n*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_tblastx_report():
    os.remove(tblastx_report)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)
