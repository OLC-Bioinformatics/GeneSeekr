#!/usr/bin/env python
from gzip import GzipFile
from StringIO import StringIO
from collections import defaultdict
import urllib
__author__ = 'mike knowles'
___doc___ = """Download CARD database and modify for use in ARMI"""

def dlunzip(url):
    compressed = StringIO()
    compressed.write(urllib.urlopen(url).read())
    #
    # Set the file's current position to the beginning
    # of the file so that gzip.GzipFile can read
    # its contents from the top.
    #
    compressed.seek(0)
    decompressed = GzipFile(fileobj=compressed, mode='rb')
    return decompressed


def makedb(path):
    genelist = defaultdict(list)
    with dlunzip('http://arpcard.mcmaster.ca/blast/db/AROtags.txt.gz') as arotags:
        for line in arotags:
            acc, name, aro = line.rstrip().split("\t")
            gene = ".".join(acc.split(".")[:-1])
            genelist[gene].append(aro[4:])
    with dlunzip('http://arpcard.mcmaster.ca/blast/db/nucleotide/ARmeta-genes.fa.gz') as g, open(path, 'w') as handle:
        for line in g:
            if ">" == line[0]:
                try:
                    aro = "".join(genelist[line.split(" ")[0][1:]])
                except KeyError:
                    aro = line.split("ARO:")[-1][:7]
                    if aro == "1000001":
                        aro = line.split("ARO:")[1][:7]
                handle.write(">ARO:%s %s" % (aro, line[1:]))
            else:
                handle.write(line)
    print 'installing CARD database into build directory'


def updatearo(path):
    import pickle
    aro = ""
    cardict = {}
    with dlunzip('http://arpcard.mcmaster.ca/obo-download/aro.obo.gz') as obo:
        for line in obo:
            if "id: ARO" == line[:7]:
                aro = line[8:].rstrip()
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
    with open(path, 'w') as handle:
        pickle.dump(cardict, handle)

if __name__ == '__main__':
    import os
    # makedb(os.path.join(os.path.split(__file__)[0], 'genes.dat'))
    updatearo(os.path.join(os.path.split(__file__)[0], 'aro.dat'))
