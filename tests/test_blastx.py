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
import blastx

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = MetadataObject()
    datapath = os.path.join(test_path, 'testdata')
    v.sequencepath = os.path.join(datapath, 'sequences')
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


blastx_method = method_init(variables(), 'geneseekr', 'blastx', True, True, blastx.BLASTx)


def test_parser():
    assert os.path.basename(blastx_method.targets[0]) == 'amr.tfa'


def test_combined_files():
    assert os.path.isfile(blastx_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastx_method.strains[0])


def test_strain():
    assert os.path.basename(blastx_method.strains[0]) == '2018-SEQ-0552.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(blastx_method.combinedtargets,
                          dbtype='prot')
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.psq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(blastx_method.metadata,
                                 blastx_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'card_aa'


def test_targetfiles():
    assert targetfiles[0] == blastx_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['yojI']


def test_blastx():
    blastx_method.metadata = geneseekr.run_blast(blastx_method.metadata,
                                                 blastx_method.analysistype,
                                                 blastx_method.program,
                                                 blastx_method.outfmt,
                                                 evalue=blastx_method.evalue,
                                                 num_threads=blastx_method.cpus)


def test_blastx_report(variables):
    global blastx_report
    blastx_report = os.path.join(variables.reportpath, '2018-SEQ-0552_blastx.csv')
    assert os.path.isfile(blastx_report)


def test_blastx_results():
    with open(blastx_report) as blast_results:
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    blastx_method.metadata = geneseekr.unique_parse_blast(blastx_method.metadata,
                                                          blastx_method.analysistype,
                                                          blastx_method.fieldnames,
                                                          blastx_method.cutoff,
                                                          blastx_method.program)
    for sample in blastx_method.metadata:
        assert sample.geneseekr.queryranges['Contig_54_76.3617'] == [[29664, 31283], [11054, 11845]]


def test_filter():
    blastx_method.metadata = geneseekr.filter_unique(blastx_method.metadata,
                                                     blastx_method.analysistype)
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    blastx_method.metadata = geneseekr.dict_initialise(blastx_method.metadata,
                                                       blastx_method.analysistype)
    for sample in blastx_method.metadata:
        assert type(sample.geneseekr.protseq) is dict


def test_report_creation():
    blastx_method.metadata = geneseekr.reporter(blastx_method.metadata,
                                                blastx_method.analysistype,
                                                blastx_method.reportpath,
                                                blastx_method.align,
                                                blastx_method.targetfiles,
                                                blastx_method.records,
                                                blastx_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(blastx_method.reportpath, 'geneseekr_blastx.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_parse_results():
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastresults['OXA_12'] == 94.19


def test_aaseq():
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastlist[0]['query_sequence'][:5] == 'MELLS' or \
               sample.geneseekr.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_combined_targets_clean():
    os.remove(blastx_method.combinedtargets)


def test_makeblastdb_clean(variables):
    databasefiles = glob(os.path.join(variables.targetpath, 'combinedtargets.p*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_blastx_report():
    os.remove(blastx_report)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)
