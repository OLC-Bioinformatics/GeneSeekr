#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.geneseekr import GeneSeekr
from geneseekr.blast import BLAST
import multiprocessing
from glob import glob
from time import time
import os

test_path = os.path.abspath(os.path.dirname(__file__))

__author__ = 'adamkoziol'


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


def method_init(analysistype, program, align, unique):
    global var
    var = variables()
    var.analysistype = analysistype
    var.program = program
    var.align = align
    var.unique = unique
    method = BLAST(var)
    return method


tblastn_method = method_init(analysistype='resfinder',
                             program='tblastn',
                             align=True,
                             unique=True)


def test_parser():
    assert os.path.basename(tblastn_method.targets[0]) == 'beta-lactam.tfa'


def test_combined_files():
    assert os.path.isfile(tblastn_method.combinedtargets)


def test_strains():
    assert os.path.isfile(tblastn_method.strains[0])


def test_strain():
    assert os.path.basename(tblastn_method.strains[0]) == 'amr_test.fasta'


def test_makeblastdb():
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(fasta=tblastn_method.combinedtargets,
                          program=tblastn_method.program)
    assert os.path.isfile(os.path.join(var.targetpath, 'combinedtargets.nsq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(metadata=tblastn_method.metadata,
                                 analysistype=tblastn_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'resfinder'


def test_targetfiles():
    assert targetfiles[0] == tblastn_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['ampH_2_HQ586946']


def test_tblastn():
    global tblastn_report
    tblastn_method.metadata = geneseekr.run_blast(metadata=tblastn_method.metadata,
                                                  analysistype=tblastn_method.analysistype,
                                                  program=tblastn_method.program,
                                                  outfmt=tblastn_method.outfmt,
                                                  evalue=tblastn_method.evalue,
                                                  num_threads=tblastn_method.cpus)
    tblastn_report = os.path.join(var.reportpath, 'amr_test_tblastn_resfinder.tsv')
    assert os.path.isfile(tblastn_report)


def test_enhance_report_parsing():
    geneseekr.parseable_blast_outputs(metadata=tblastn_method.metadata,
                                      analysistype=tblastn_method.analysistype,
                                      fieldnames=tblastn_method.fieldnames,
                                      program=tblastn_method.program)
    header = open(tblastn_report).readline()
    assert header.split('\t')[0] == 'query_id'


def test_tblastn_results():
    with open(tblastn_report) as blast_results:
        next(blast_results)
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    tblastn_method.metadata = geneseekr.unique_parse_blast(metadata=tblastn_method.metadata,
                                                           analysistype=tblastn_method.analysistype,
                                                           fieldnames=tblastn_method.fieldnames,
                                                           cutoff=tblastn_method.cutoff,
                                                           program=tblastn_method.program)
    for sample in tblastn_method.metadata:
        assert sample.resfinder.queryranges['contig2'] == [[1, 264]]


def test_filter():
    tblastn_method.metadata = geneseekr.filter_unique(metadata=tblastn_method.metadata,
                                                      analysistype=tblastn_method.analysistype)
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    tblastn_method.metadata = geneseekr.dict_initialise(metadata=tblastn_method.metadata,
                                                        analysistype=tblastn_method.analysistype)
    for sample in tblastn_method.metadata:
        assert type(sample.resfinder.protseq) is dict


def test_report_creation():
    tblastn_method.metadata = geneseekr.resfinder_reporter(metadata=tblastn_method.metadata,
                                                           analysistype=tblastn_method.analysistype,
                                                           reportpath=tblastn_method.reportpath,
                                                           align=tblastn_method.align,
                                                           program=tblastn_method.program,
                                                           targetpath=tblastn_method.targetpath,
                                                           cutoff=tblastn_method.cutoff)


def test_report_existance():
    global geneseekr_report
    geneseekr_report = os.path.join(tblastn_method.reportpath, 'resfinder_tblastn.xlsx')
    assert os.path.isfile(geneseekr_report)


def test_report_row():
    for sample in tblastn_method.metadata:
        assert sorted(sample.resfinder.sampledata)[0][0] == 'blaOXA'


def test_parse_results():
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastresults['blaOXA_427_1_KX827604'] == 94.34


def test_aaseq():
    for sample in tblastn_method.metadata:
        assert sample.resfinder.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_fasta_create():
    global fasta_file
    geneseekr.export_fasta(metadata=tblastn_method.metadata,
                           analysistype=tblastn_method.analysistype,
                           reportpath=tblastn_method.reportpath,
                           cutoff=tblastn_method.cutoff,
                           program=tblastn_method.program)
    fasta_file = os.path.join(var.reportpath, 'amr_test_resfinder.fasta')
    assert os.path.isfile(fasta_file)
    header = open(fasta_file, 'r').readline().rstrip()
    assert header == '>amr_test_blaOXA_427_1_KX827604'


def test_combined_targets_clean():
    os.remove(tblastn_method.combinedtargets)


def test_makeblastdb_clean():
    databasefiles = glob(os.path.join(var.targetpath, 'combinedtargets.n*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_tblastn_report():
    os.remove(tblastn_report)


def test_remove_fasta_file():
    os.remove(fasta_file)


def test_remove_geneseekr_report():
    os.remove(geneseekr_report)


def test_remove_report_path():
    os.rmdir(tblastn_method.reportpath)
