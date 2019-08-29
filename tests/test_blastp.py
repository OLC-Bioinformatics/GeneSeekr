#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import MetadataObject
from genemethods.geneseekr.geneseekr import GeneSeekr
from genemethods.geneseekr.blast import BLAST
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


def method_init(analysistype, program, align, unique):
    global var
    var = variables()
    var.analysistype = analysistype
    var.program = program
    var.align = align
    var.unique = unique
    method = BLAST(var)
    return method


blastp_method = method_init(analysistype='geneseekr',
                            program='blastp',
                            align=True,
                            unique=True)


def test_parser():
    assert os.path.basename(blastp_method.targets[0]) == 'amr.tfa'


def test_combined_files():
    assert os.path.isfile(blastp_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastp_method.strains[0])


def test_strain():
    assert os.path.basename(blastp_method.strains[0]) == 'amr_test.fasta'


def test_makeblastdb():
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(fasta=blastp_method.combinedtargets,
                          program=blastp_method.program)
    assert os.path.isfile(os.path.join(var.targetpath, 'combinedtargets.psq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(metadata=blastp_method.metadata,
                                 analysistype=blastp_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'card_aa'


def test_targetfiles():
    assert targetfiles[0] == blastp_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['yojI']


def test_blastp():
    global blastp_report
    blastp_method.metadata = geneseekr.run_blast(metadata=blastp_method.metadata,
                                                 analysistype=blastp_method.analysistype,
                                                 program=blastp_method.program,
                                                 outfmt=blastp_method.outfmt,
                                                 evalue=blastp_method.evalue,
                                                 num_threads=blastp_method.cpus)
    blastp_report = os.path.join(var.reportpath, 'amr_test_blastp_geneseekr.tsv')
    assert os.path.isfile(blastp_report)


def test_enhance_report_parsing():
    geneseekr.parseable_blast_outputs(metadata=blastp_method.metadata,
                                      analysistype=blastp_method.analysistype,
                                      fieldnames=blastp_method.fieldnames,
                                      program=blastp_method.program)
    header = open(blastp_report).readline()
    assert header.split('\t')[0] == 'query_id'


def test_blastp_results():
    with open(blastp_report) as blast_results:
        next(blast_results)
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    blastp_method.metadata = geneseekr.unique_parse_blast(metadata=blastp_method.metadata,
                                                          analysistype=blastp_method.analysistype,
                                                          fieldnames=blastp_method.fieldnames,
                                                          cutoff=blastp_method.cutoff,
                                                          program=blastp_method.program)
    for sample in blastp_method.metadata:
        assert sample.geneseekr.queryranges['contig1'] == [[1, 547]]


def test_filter():
    blastp_method.metadata = geneseekr.filter_unique(metadata=blastp_method.metadata,
                                                     analysistype=blastp_method.analysistype)
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    blastp_method.metadata = geneseekr.dict_initialise(metadata=blastp_method.metadata,
                                                       analysistype=blastp_method.analysistype)
    for sample in blastp_method.metadata:
        assert type(sample.geneseekr.protseq) is dict


def test_report_creation():
    blastp_method.metadata = geneseekr.resfinder_reporter(metadata=blastp_method.metadata,
                                                          analysistype=blastp_method.analysistype,
                                                          reportpath=blastp_method.reportpath,
                                                          align=blastp_method.align,
                                                          program=blastp_method.program,
                                                          targetpath=blastp_method.targetpath,
                                                          cutoff=blastp_method.cutoff)


def test_report_csv():
    global geneseekr_csv
    geneseekr_csv = os.path.join(blastp_method.reportpath, 'amr_test_blastp_geneseekr.tsv')
    assert os.path.isfile(geneseekr_csv)


def test_report_xls():
    global geneseekr_xls
    geneseekr_xls = os.path.join(blastp_method.reportpath, 'geneseekr_blastp.xlsx')
    assert os.path.isfile(geneseekr_xls)


def test_parse_results():
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastresults['OXA_12'] == 91.86


def test_aaseq():
    for sample in blastp_method.metadata:
        assert sample.geneseekr.blastlist[0]['query_sequence'][:4] == 'MELL' or \
               sample.geneseekr.blastlist[0]['query_sequence'][:4] == 'MSRI'


def test_fasta_create():
    global fasta_file
    geneseekr.export_fasta(metadata=blastp_method.metadata,
                           analysistype=blastp_method.analysistype,
                           reportpath=blastp_method.reportpath,
                           cutoff=blastp_method.cutoff,
                           program=blastp_method.program)
    fasta_file = os.path.join(var.reportpath, 'amr_test_geneseekr.fasta')
    assert os.path.isfile(fasta_file)
    header = open(fasta_file, 'r').readline().rstrip()
    assert header == '>amr_test_OXA_12'


def test_combined_targets_clean():
    os.remove(blastp_method.combinedtargets)


def test_makeblastdb_clean():
    databasefiles = glob(os.path.join(var.targetpath, 'combinedtargets.p*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_blastp_report():
    os.remove(blastp_report)


def test_remove_fasta_file():
    os.remove(fasta_file)


def test_remove_geneseekr_xls():
    os.remove(geneseekr_xls)


def test_remove_report_path():
    os.rmdir(blastp_method.reportpath)
