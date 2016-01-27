#!/usr/bin/env python
from ARMI_Lt import ARMISeekr
from Bio.Blast.Applications import NcbiblastnCommandline
__author__ = 'mike knowles'


class Jackson(ARMISeekr):

    def _blast(self, (fasta, db)):
        blastn = NcbiblastnCommandline(query=fasta,
                                       db=db,
                                       evalue=10,
                                       outfmt="'6 sseqid nident slen'",
                                       perc_identity=self.cutoff,
                                       num_descriptions=10000,
                                       num_alignments=10000)
        stdout, stderr = blastn()
        if stdout != '':
            return [[[fasta], [aln[0],
                     (lambda x: "{:.2f}%".format(x*100) if x < 1.0 else '+')(abs(float(aln[1]) / float(aln[2])))]]
                    for aln in [hsp.split('\t')
                    for hsp in stdout.rstrip().split("\n")]
                    if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]

def resfinder(ardb, plus, out):
    from collections import defaultdict
    import os
    import time
    import json
    assert isinstance(out, str), u'Output location is not a string {0!r:s}'.format(out)
    assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
    print "[{}] Writing CSV and JSON to output directory".format(time.strftime("%H:%M:%S"))
    resistlist = set(r for gene in ardb if 'resist' in ardb[gene] for r in ardb[gene]['resist'])
    resistlist.add('tellurite')
    resistlist = sorted(resistlist)
    rowcount, row = 0, 'Strain,'
    row += ', '.join(resistlist)
    anti = dict((genome, dict()) for genome in plus)
    for genomerow in sorted(plus):
        lowlst = [x.lower() for x in plus[genomerow]]
        row += '\n{},'.format(os.path.split(os.path.splitext(genomerow)[0])[1].replace('_filteredAssembled', ""))
        for gene in plus[genomerow]:
            req = True
            if gene in ('TehA', 'TehB'):
                anti[genomerow]['tellurite'] = [{gene: plus[genomerow][gene]}]
            try:
                if 'requires' in ardb[gene.lower()]:
                    test = len(ardb[gene.lower()]['requires'])
                    if len([x for x in ardb[gene.lower()]['requires'] if x in lowlst]) >= test:
                        req = True
                    else:
                        req = False
                if 'resist' in ardb[gene.lower()] and req:
                    for resist in ardb[gene.lower()]['resist']:
                        if resist not in anti[genomerow]:
                            anti[genomerow][resist] = []
                        anti[genomerow][resist].append({gene: plus[genomerow][gene]})
            except KeyError:
                print gene
        for resist in resistlist:
            if resist in anti[genomerow]:
                row += '\"{0}\",'.format(','.join(sorted(set(y for x in anti[genomerow][resist] for y in x))))
            else:
                row += '-,'
            # row += ',' + (lambda x, y: ' '.join(map(str, x[y])) if y in x else 'N')(plus[genomerow], genename)
            # Add the allele numbers to the row for the appropriate gene, otherwise return N
    with open("%s/resfinder_results_%s.csv" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(row)
    with open("%s/resfinder_results_%s.json" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as j:
        json.dump(anti, j, sort_keys=True,indent=4, separators=(',', ': '))
