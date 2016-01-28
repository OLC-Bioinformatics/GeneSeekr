#!/usr/bin/env python
from gzip import GzipFile
from StringIO import StringIO
from collections import defaultdict
import urllib
__author__ = 'mike knowles'


def dlunzip(url):
    source = urllib.urlopen(url).read()
    compressedFile = StringIO()
    compressedFile.write(source)
    #
    # Set the file's current position to the beginning
    # of the file so that gzip.GzipFile can read
    # its contents from the top.
    #
    compressedFile.seek(0)
    decompressedFile = GzipFile(fileobj=compressedFile, mode='rb')
    return decompressedFile

def makedb(path):
    print os.getcwd()
    genelist = defaultdict(list)
    for line in dlunzip('http://arpcard.mcmaster.ca/blast/db/AROtags.txt.gz'):
        acc, name, aro = line.rstrip().split("\t")
        gene = ".".join(acc.split(".")[:-1])
        genelist[gene].append(aro[4:])
    handle = open(path, 'w')
    for line in dlunzip('http://arpcard.mcmaster.ca/blast/db/nucleotide/ARmeta-genes.fa.gz'):
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
    handle.close()
    print 'installing card database'
if __name__ == '__main__':
    import os
    makedb(os.path.join(os.path.split(__file__)[0], 'ARMI-genes.dat'))