#!/usr/bin/env python
import os
import pysamstats
from ARMI_Lt import ARMISeekr, KeyboardInterruptError, lcs
from subprocess import Popen, PIPE, STDOUT
from bowtie import Bowtie2CommandLine, Bowtie2BuildCommandLine
from Bio.Sequencing.Applications import SamtoolsViewCommandline, SamtoolsSortCommandline, SamtoolsIndexCommandline

__author__ = 'mike knowles'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        import errno
        if exception.errno != errno.EEXIST:
            raise


def get_version(exe):
    """
    :param exe: :type list required
    """
    assert isinstance(exe, list)
    return Popen(exe, stdout=PIPE, stderr=STDOUT).stdout.read()


class RawARMI(ARMISeekr):

    def __init__(self, subject, query, **kwargs):

        super(RawARMI, self).__init__(subject, query, aligner="bowtie2", **kwargs)
        self.db = subject  # remove the file extension for easier globing

    def makeblastdb(self, (fasta, db)):
        try:
            if not os.path.isfile('{}.1.bt2'.format(db)) or self.recreate:  # add nhr for searching
                assert os.path.isfile(fasta)  # check that the fasta has been specified properly
                Bowtie2BuildCommandLine(f=True,
                                        threads=self.threads,
                                        bt2=db,
                                        reference=fasta)()
            return 0
        except KeyboardInterrupt:
            raise KeyboardInterruptError()

    def _bowtie(self, raw, db):
        version = Popen(['samtools', '--version'], stdout=PIPE, stderr=STDOUT).stdout.read().split('\n')[0].split()[1]
        raw = map(os.path.abspath, raw)
        if len(raw) == 2:
            name = lcs(*raw)
            indict = dict(("m" + str(x), fastq) for x, fastq in enumerate(raw, 1))
        else:
            indict = dict(("U", ",".join(raw)))
            name = os.path.splitext(raw)[0]
        # SAMtools sort v1.3 has different run parameters
        workingdir = name + "tmp"
        make_path(workingdir)
        name += ".sorted.bam"
        if version < "1.3":
            samsort = SamtoolsSortCommandline(input_bam="-", out_prefix=name)
        else:
            samsort = SamtoolsSortCommandline(input_bam=name, o=True, out_prefix="-")
        indict.update(dict(samtools=[SamtoolsViewCommandline(b=True, S=True, input_file="-"), samsort]))
        if not os.path.isfile(name):
            Bowtie2CommandLine(bt2=os.path.splitext(os.path.abspath(db))[0],
                               threads=self.threads,
                               very_sensitive_local=True,
                               a=True,
                               **indict)(cwd=workingdir)
        if not os.path.isfile(name + ".bai"):
            SamtoolsIndexCommandline(input_bam=name)(cwd=workingdir)
        os.rmdir(workingdir)
        genes = dict()
        for rec in pysamstats.stat_baseq_ext(alignmentfile=name, fafile=db):
            # Values of interest can be retrieved using the appropriate keys
            # Simple filtering statement: if the number of matches at a particular position in the reference sequence is
            # greater than the number of mismatches, and the total depth is 5 or more, add the position of the results
            if rec['chrom'] not in genes:
                genes[rec['chrom']] = dict(identity=1.0, depth=float(rec['reads_all']))
            else:
                genes[rec['chrom']]['identity'] += 1.0
                genes[rec['chrom']]['depth'] += float(rec['reads_all'])
        return genes

    def _blast(self, (raw, db)):
        genes = self._bowtie(raw, db)
        for gene in genes:
            db, accn, span, aro, name = gene.split('|')

            start, end = map(int, span.split('-'))
            length = end - start + 1
            outof = "{0:d}/{1:d}".format(*map(int, [genes[gene]['identity'], length]))
            genes[gene]['identity'] /= float(length) / 100
            genes[gene]['depth'] /= float(length)
            # Add cutoff
            if genes[gene]['identity'] >= self.cutoff and genes[gene]['depth'] > 4.0:
                yield [[aro[4:]],
                       [genes[gene]['identity'],
                        'average depth:' + str(genes[gene]['depth']),
                        db, accn, span, outof, name]]


if __name__ == '__main__':
    from argparse import ArgumentParser
    import json
    parent = ArgumentParser()
    parent.add_argument('--subject', default=['data/genes.dat'], help='Specify input fasta folder')
    # parent.add_argument('--query', default=[os.getcwd()], help='Specify output folder for csv and json')
    parent.add_argument('-t', '--threads', type=int, default=4, help='Specify number of threads')
    parent.add_argument('--recreate', action='store_true', help='recreate alignment databases')
    from glob import iglob
    # use the lcs to determine the fastq pairs
    quandry = [y for y in iglob('/Users/mike/Documents/armi/EHEC-SalmRT-LmonoHist-26528831/2015-SEQ-1172-30853222/Data/'
                                'Intensities/BaseCalls/*.gz')]
    test = RawARMI(query=[quandry], **vars(parent.parse_args()))
    with open('/Users/mike/Documents/armi/output.json', 'w') as out:
        json.dump(test.mpblast(cutoff=95), out, indent=4, sort_keys=True)
