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
import blastn

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
def method_init(variables, analysistype, program, dest):
    global method
    variables.analysistype = analysistype
    variables.program = program
    method = dest(variables)
    return method


blastn_method = method_init(variables(), 'geneseekr', 'blastn', blastn.BLASTn)


def test_parser():
    assert os.path.basename(blastn_method.targets[0]) == 'aminoglycoside.tfa'


def test_combined_files():
    assert os.path.isfile(blastn_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastn_method.strains[0])


def test_strain():
    assert os.path.basename(blastn_method.strains[0]) == '2018-SEQ-0552.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(blastn_method.combinedtargets)
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.nsq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(blastn_method.metadata,
                                 blastn_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'resfinder'


def test_targetfiles():
    assert targetfiles[0] == blastn_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['blaOXA_235_1_JQ820240']


def test_blastn():
    blastn_method.metadata = geneseekr.run_blast(blastn_method.metadata,
                                                 blastn_method.analysistype,
                                                 blastn_method.program,
                                                 blastn_method.outfmt,
                                                 evalue=blastn_method.evalue,
                                                 num_threads=blastn_method.cpus)


def test_blastn_report(variables):
    global blastn_report
    blastn_report = os.path.join(variables.reportpath, '2018-SEQ-0552_blastn.csv')
    assert os.path.isfile(blastn_report)


def test_blastn_results():
    with open(blastn_report) as blast_results:
        data = blast_results.readline()
        results = data.split('\t')
        assert results[2] == '179'


def test_blast_parse():
    blastn_method.metadata = geneseekr.parse_blastn(blastn_method.metadata,
                                                    blastn_method.analysistype,
                                                    blastn_method.fieldnames,
                                                    blastn_method.cutoff)


def test_report_creation():
    blastn_method.metadata = geneseekr.reporter(blastn_method.metadata,
                                                blastn_method.analysistype,
                                                blastn_method.reportpath,
                                                blastn_method.align,
                                                blastn_method.targetfiles,
                                                blastn_method.records,
                                                blastn_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(blastn_method.reportpath, 'geneseekr.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_parse_results():
    for sample in blastn_method.metadata:
        assert sample.geneseekr.blastresults['cphA6_1_AY227052'] == 88.1


def test_combined_targets_clean():
    os.remove(blastn_method.combinedtargets)


def test_makeblastdb_clean(variables):
    databasefiles = glob(os.path.join(variables.targetpath, 'combinedtargets.n*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_blastn_report():
    os.remove(blastn_report)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)
