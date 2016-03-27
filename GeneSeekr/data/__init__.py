#!/usr/bin/env python
import tarfile
import urllib
from StringIO import StringIO
from collections import defaultdict
__author__ = 'mike knowles'
___doc___ = """Download CARD database and modify for use in ARMI"""

class Build:

    def __init__(self):
        self.compressed = StringIO()
        self.compressed.write(urllib.urlopen('https://card.mcmaster.ca/download/0/broadsteet-v1.0.3.tar.gz').read())
        #
        # Set the file's current position to the beginning
        # of the file so that gzip.GzipFile can read
        # its contents from the top.
        #
        self.compressed.seek(0)
        self.tar = tarfile.open(fileobj=self.compressed)
        self.fasta = self.tar.getmember('./nucleotide_fasta[protein homolog model].fasta')
        self.aro = self.tar.getmember('./aro.obo')

    def updatearo(self, path):
        import pickle
        aro = ""
        cardict = {}
        for line in self.tar.extractfile(self.aro):
            if "id: ARO" == line[:7]:
                aro = line[8:].rstrip()
                if aro not in cardict:
                    cardict[aro] = defaultdict(list)
            if "name: " == line[:6]:
                cardict[aro]["name"] = line[6:].rstrip()
            if "is_a: ARO:" in line:
                cardict[aro]["isa"].append(line[10:18].rstrip())
                cardict[aro]["function"] = line[20:].rstrip()
            if "relationship: confers_resistance_to_drug" in line:
                cardict[aro]["resist"].append(line[55:].rstrip())
            elif "relationship: confers_resistance_to " in line:
                cardict[aro]["resist"].append(line[50:].rstrip())
            if "relationship: part_of" in line:
                cardict[aro]["complex"].append(line[26:33].rstrip())
                if line[26:33].rstrip() not in cardict:
                    cardict[line[26:33].rstrip()] = defaultdict(list)
                cardict[line[26:33].rstrip()]["member"].append(aro)
                cardict[line[26:33].rstrip()]["member"].sort()
            if "relationship: targeted_by_drug" in line:
                cardict[aro]["sensitivity"].append(line[45:].rstrip())
        print "updating CARD gene onotologies"
        with open(path, 'w') as handle:
            pickle.dump(cardict, handle)

    def makedb(self, path):
        print 'installing CARD database into build directory'
        self.fasta.name = path
        self.tar.extract(self.fasta)

    def __del__(self):
        self.tar.close()
        self.compressed.close()
        del self


if __name__ == '__main__':
    card = Build()
    import os
    card.makedb(os.path.join(os.path.split(__file__)[0], 'genes.dat'))
    card.updatearo(os.path.join(os.path.split(__file__)[0], 'aro.dat'))
    del card
