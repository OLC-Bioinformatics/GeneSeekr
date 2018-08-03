#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.geneseekr import BLAST, GeneSeekr
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
def method_init(variables, analysistype, program, align, unique):
    global method
    variables.analysistype = analysistype
    variables.program = program
    variables.align = align
    variables.unique = unique
    method = BLAST(variables)
    return method


blastn_method = method_init(variables(), 'resfinder', 'blastn', True, True)


def test_parser():
    assert os.path.basename(blastn_method.targets[0]) == 'beta-lactam.tfa'


def test_combined_files():
    assert os.path.isfile(blastn_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastn_method.strains[0])


def test_strain():
    assert os.path.basename(blastn_method.strains[0]) == '2018-SEQ-0552.fasta'


def test_makeblastdb(variables):
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(blastn_method.combinedtargets,
                          blastn_method.program)
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
    assert records[targetfiles[0]]['blaOXA_427_1_KX827604']


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
        assert int(results[2]) >= 179


def test_blast_parse():
    blastn_method.metadata = geneseekr.unique_parse_blast(blastn_method.metadata,
                                                          blastn_method.analysistype,
                                                          blastn_method.fieldnames,
                                                          blastn_method.cutoff,
                                                          blastn_method.program)
    for sample in blastn_method.metadata:
        assert sample.resfinder.queryranges['Contig_54_76.3617'] == [[11054, 11848]]


def test_filter():
    blastn_method.metadata = geneseekr.filter_unique(blastn_method.metadata,
                                                     blastn_method.analysistype)
    for sample in blastn_method.metadata:
        assert sample.resfinder.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    blastn_method.metadata = geneseekr.dict_initialise(blastn_method.metadata,
                                                       blastn_method.analysistype)
    for sample in blastn_method.metadata:
        assert type(sample.resfinder.protseq) is dict

def test_report_creation():
    blastn_method.metadata = geneseekr.resfinder_reporter(blastn_method.metadata,
                                                          blastn_method.analysistype,
                                                          targetfolders,
                                                          blastn_method.reportpath,
                                                          blastn_method.align,
                                                          blastn_method.targetfiles,
                                                          blastn_method.records,
                                                          blastn_method.program)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(blastn_method.reportpath, 'resfinder_blastn.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_report_row():
    for sample in blastn_method.metadata:
        assert sorted(sample.resfinder.sampledata)[0] == \
               ['blaOXA', '1', 'Beta-Lactamase', 86.16, 100.0, 'Contig_54_76.3617', '11054...11848', '-']


def test_parse_results():
    for sample in blastn_method.metadata:
        assert sample.resfinder.blastresults['blaOXA_427_1_KX827604'] == 86.16


def test_aaseq():
    for sample in blastn_method.metadata:
        assert sample.resfinder.protseq['blaOXA_427_1_KX827604'][:5] == 'MSRIL'


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


def test_remove_report_path():
    os.rmdir(blastn_method.reportpath)
