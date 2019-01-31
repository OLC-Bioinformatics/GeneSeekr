#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.geneseekr import GeneSeekr
from geneseekr.blast import BLAST
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


tblastx_method = method_init(variables(), 'resfinder', 'tblastx', True, True)


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
    geneseekr.makeblastdb(fasta=tblastx_method.combinedtargets,
                          program=tblastx_method.program)
    assert os.path.isfile(os.path.join(variables.targetpath, 'combinedtargets.nsq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(metadata=tblastx_method.metadata,
                                 analysistype=tblastx_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'resfinder'


def test_targetfiles():
    assert targetfiles[0] == tblastx_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['ampH_2_HQ586946']


def test_tblastx(variables):
    global tblastx_report
    tblastx_method.metadata = geneseekr.run_blast(metadata=tblastx_method.metadata,
                                                  analysistype=tblastx_method.analysistype,
                                                  program=tblastx_method.program,
                                                  outfmt=tblastx_method.outfmt,
                                                  evalue=tblastx_method.evalue,
                                                  num_threads=tblastx_method.cpus)
    tblastx_report = os.path.join(variables.reportpath, '2018-SEQ-0552_tblastx_resfinder.tsv')
    assert os.path.isfile(tblastx_report)



def test_enhance_report_parsing():
    geneseekr.parseable_blast_outputs(metadata=tblastx_method.metadata,
                                      analysistype=tblastx_method.analysistype,
                                      fieldnames=tblastx_method.fieldnames,
                                      program=tblastx_method.program)
    header = open(tblastx_report).readline()
    assert header.split('\t')[0] == 'query_id'


def test_tblastx_results():
    with open(tblastx_report) as blast_results:
        next(blast_results)
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    tblastx_method.metadata = geneseekr.unique_parse_blast(metadata=tblastx_method.metadata,
                                                           analysistype=tblastx_method.analysistype,
                                                           fieldnames=tblastx_method.fieldnames,
                                                           cutoff=tblastx_method.cutoff,
                                                           program=tblastx_method.program)
    for sample in tblastx_method.metadata:
        assert sample.resfinder.queryranges['Contig_54_76.3617'] == [[11054, 11848]]


def test_filter():
    tblastx_method.metadata = geneseekr.filter_unique(metadata=tblastx_method.metadata,
                                                      analysistype=tblastx_method.analysistype)
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    tblastx_method.metadata = geneseekr.dict_initialise(metadata=tblastx_method.metadata,
                                                        analysistype=tblastx_method.analysistype)
    for sample in tblastx_method.metadata:
        assert type(sample.resfinder.protseq) is dict


def test_report_creation():
    tblastx_method.metadata = geneseekr.resfinder_reporter(metadata=tblastx_method.metadata,
                                                           analysistype=tblastx_method.analysistype,
                                                           reportpath=tblastx_method.reportpath,
                                                           align=tblastx_method.align,
                                                           program=tblastx_method.program,
                                                           targetpath=tblastx_method.targetpath,
                                                           cutoff=tblastx_method.cutoff)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(tblastx_method.reportpath, 'resfinder_tblastx.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_report_row():
    for sample in tblastx_method.metadata:
        assert sorted(sample.resfinder.sampledata)[0][0] == 'ampH'


def test_parse_results():
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastresults['ampH_2_HQ586946'] == 93.96


def test_aaseq():
    for sample in tblastx_method.metadata:
        assert sample.resfinder.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_fasta_create(variables):
    global fasta_file
    geneseekr.export_fasta(metadata=tblastx_method.metadata,
                           analysistype=tblastx_method.analysistype,
                           reportpath=tblastx_method.reportpath,
                           cutoff=tblastx_method.cutoff,
                           program=tblastx_method.program)
    fasta_file = os.path.join(variables.reportpath, '2018-SEQ-0552_resfinder.fasta')
    assert os.path.isfile(fasta_file)
    header = open(fasta_file, 'r').readline().rstrip()
    assert header == '>2018-SEQ-0552_ampH_2_HQ586946'


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


def test_remove_fasta_file():
    os.remove(fasta_file)

def test_remove_report_path():
    os.rmdir(tblastx_method.reportpath)
