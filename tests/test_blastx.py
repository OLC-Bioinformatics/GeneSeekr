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


def method_init(analysistype, program, align, unique):
    global var
    var = variables()
    var.analysistype = analysistype
    var.program = program
    var.align = align
    var.unique = unique
    method = BLAST(var)
    return method


blastx_method = method_init(analysistype='geneseekr',
                            program='blastx',
                            align=True,
                            unique=True)


def test_parser():
    assert os.path.basename(blastx_method.targets[0]) == 'amr.tfa'


def test_combined_files():
    assert os.path.isfile(blastx_method.combinedtargets)


def test_strains():
    assert os.path.isfile(blastx_method.strains[0])


def test_strain():
    assert os.path.basename(blastx_method.strains[0]) == '2018-SEQ-0552.fasta'


def test_makeblastdb():
    global geneseekr
    geneseekr = GeneSeekr()
    geneseekr.makeblastdb(fasta=blastx_method.combinedtargets,
                          program=blastx_method.program)
    assert os.path.isfile(os.path.join(var.targetpath, 'combinedtargets.psq'))


def test_variable_populate():
    global targetfolders
    global targetfiles
    global records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(metadata=blastx_method.metadata,
                                 analysistype=blastx_method.analysistype)


def test_targetfolders():
    assert os.path.basename(list(targetfolders)[0]) == 'card_aa'


def test_targetfiles():
    assert targetfiles[0] == blastx_method.combinedtargets


def test_records():
    assert records[targetfiles[0]]['yojI']


def test_blastx():
    global blastx_report
    blastx_method.metadata = geneseekr.run_blast(metadata=blastx_method.metadata,
                                                 analysistype=blastx_method.analysistype,
                                                 program=blastx_method.program,
                                                 outfmt=blastx_method.outfmt,
                                                 evalue=blastx_method.evalue,
                                                 num_threads=blastx_method.cpus)
    blastx_report = os.path.join(var.reportpath, '2018-SEQ-0552_blastx_geneseekr.tsv')
    assert os.path.isfile(blastx_report)


def test_enhance_report_parsing():
    geneseekr.parseable_blast_outputs(metadata=blastx_method.metadata,
                                      analysistype=blastx_method.analysistype,
                                      fieldnames=blastx_method.fieldnames,
                                      program=blastx_method.program)
    header = open(blastx_report).readline()
    assert header.split('\t')[0] == 'query_id'


def test_blastx_results():
    with open(blastx_report) as blast_results:
        next(blast_results)
        data = blast_results.readline()
        results = data.split('\t')
        assert int(results[2]) >= 50


def test_blast_parse():
    blastx_method.metadata = geneseekr.unique_parse_blast(metadata=blastx_method.metadata,
                                                          analysistype=blastx_method.analysistype,
                                                          fieldnames=blastx_method.fieldnames,
                                                          cutoff=blastx_method.cutoff,
                                                          program=blastx_method.program)
    for sample in blastx_method.metadata:
        assert sample.geneseekr.queryranges['Contig_54_76.3617'] == [[29664, 31283], [11054, 11845]]


def test_filter():
    blastx_method.metadata = geneseekr.filter_unique(metadata=blastx_method.metadata,
                                                     analysistype=blastx_method.analysistype)
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastlist[0]['percentidentity'] >= 70


def test_dict_create():
    blastx_method.metadata = geneseekr.dict_initialise(metadata=blastx_method.metadata,
                                                       analysistype=blastx_method.analysistype)
    for sample in blastx_method.metadata:
        assert type(sample.geneseekr.protseq) is dict


def test_target_folders():
    global targetfolders, targetfiles, records
    targetfolders, targetfiles, records = \
        geneseekr.target_folders(metadata=blastx_method.metadata,
                                 analysistype=blastx_method.analysistype)
    assert records[targetfiles[0]]['yojI']


def test_report_creation():
    blastx_method.metadata = geneseekr.reporter(metadata=blastx_method.metadata,
                                                analysistype=blastx_method.analysistype,
                                                reportpath=blastx_method.reportpath,
                                                align=blastx_method.align,
                                                records=records,
                                                program=blastx_method.program,
                                                cutoff=blastx_method.cutoff)


def test_report_csv():
    global geneseekr_csv
    geneseekr_csv = os.path.join(blastx_method.reportpath, 'geneseekr_blastx.csv')
    assert os.path.isfile(geneseekr_csv)


def test_detailed_report_csv():
    global geneseekr_detailed_csv
    geneseekr_detailed_csv = os.path.join(blastx_method.reportpath, 'geneseekr_blastx_detailed.csv')
    assert os.path.isfile(geneseekr_detailed_csv)


def test_report_xls():
    global geneseekr_xls
    geneseekr_xls = os.path.join(blastx_method.reportpath, 'geneseekr_blastx.xlsx')
    assert os.path.isfile(geneseekr_xls)


def test_parse_results():
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastresults['OXA_12'] == 91.86


def test_aaseq():
    for sample in blastx_method.metadata:
        assert sample.geneseekr.blastlist[0]['query_sequence'][:5] == 'MELLS' or \
               sample.geneseekr.blastlist[0]['query_sequence'][:5] == 'MSRIL'


def test_fasta_create():
    global fasta_file
    geneseekr.export_fasta(metadata=blastx_method.metadata,
                           analysistype=blastx_method.analysistype,
                           reportpath=blastx_method.reportpath,
                           cutoff=blastx_method.cutoff,
                           program=blastx_method.program)
    fasta_file = os.path.join(var.reportpath, '2018-SEQ-0552_geneseekr.fasta')
    assert os.path.isfile(fasta_file)
    header = open(fasta_file, 'r').readline().rstrip()
    assert header == '>2018-SEQ-0552_OXA_12'


def test_combined_targets_clean():
    os.remove(blastx_method.combinedtargets)


def test_makeblastdb_clean():
    databasefiles = glob(os.path.join(var.targetpath, 'combinedtargets.p*'))
    for dbfile in databasefiles:
        os.remove(dbfile)


def test_remove_blastx_report():
    os.remove(blastx_report)


def test_remove_geneseekr_csv():
    os.remove(geneseekr_csv)


def test_remove_fasta_file():
    os.remove(fasta_file)


def test_removed_detailed_geneseekr_csv():
    os.remove(geneseekr_detailed_csv)


def test_remove_geneseekr_xls():
    os.remove(geneseekr_xls)


def test_remove_report_path():
    os.rmdir(blastx_method.reportpath)
