#! /usr/env/python
__author__ = 'mikeknowles, akoziol'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
Revised with speed improvements
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
from cStringIO import StringIO
from glob import glob
import subprocess, os, time, sys, shlex, re, threading, json, mmap, errno
from argparse import ArgumentParser

parser = ArgumentParser(description='Performs blast analyses to determine presence of custom targets')
parser.add_argument('-p', '--path', required=False,
                    default='/home/blais/PycharmProjects/pythonGeneSeekr/',
                    help='Specify path for custom folder locations')
parser.add_argument('-c', '--cutoff', required=False, default=0.8,
                    help='The identity cutoff value for BLAST matches. Default is 0.8')
parser.add_argument('-s', '--sequencePath', required=False,
                    default='/home/blais/PycharmProjects/pythonGeneSeekr/sequences',
                    help='The location of the query sequence files')
parser.add_argument('-t', '--targetPath', required=False,
                    default='/home/blais/PycharmProjects/pythonGeneSeekr/Organism',
                    help='The location of the target files')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
# Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
path = os.path.join(args['path'], "")
cutoff = float(args['cutoff'])
sequencePath = os.path.join(args['sequencePath'], "")
targetPath = os.path.join(args['targetPath'], "")

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the count used in the dotter function
count = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def makeblastdb(dqueue):
    while True:  # while daemon
        fastapath = dqueue.get()  # grabs fastapath from dqueue
        # remove the path and the file extension for easier future globbing
        db = fastapath.split(".")[0]
        nhr = "%s.nhr" % db  # add nhr for searching
        # print nhr
        FNULL = open(os.devnull, 'w')  # define /dev/null
        if not os.path.isfile(str(nhr)):  # if check for already existing dbs
            subprocess.call(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, db)), stdout=FNULL, stderr=FNULL)
            # make blastdb
            dotter()
        dqueue.task_done()  # signals to dqueue job is done
        sys.exit()  # necessary

# Declare queues, list, and dict
dqueue = Queue()
blastqueue = Queue()
parsequeue = Queue()
testqueue = Queue()
plusqueue = Queue()
plusdict = {}
genedict = defaultdict(list)
blastpath = []
# threadlock, useful for overcoming GIL
threadlock = threading.Lock()


def makedbthreads(fastas):
    """Setup and create threads for class"""
    # Create and start threads for each fasta file in the list
    for i in range(len(fastas)):
        # Send the threads to makeblastdb
        threads = Thread(target=makeblastdb, args=(dqueue,))
        # Set the daemon to true - something to do with thread management
        threads.setDaemon(True)
        # Start the threading
        threads.start()
    for fasta in fastas:
        # Add the fasta file to the queue
        dqueue.put(fasta)
    dqueue.join()  # wait on the dqueue until everything has been processed


def xmlout(fasta, genome):
    """Parses variables from supplied tuples? dictionaries?"""
    global path
    pathString = genome.split("/")
    # Create the necessary variables from the supplied strings
    gene = fasta.split('/')[-1]  # split file from path, could use os.path.split
    genename = gene.split('.')[0]
    genomename = pathString[-1].split('.')[0]
    # Create the out variable containing the path and name of BLAST output file
    tmpPath = "%stmp" % path
    make_path(tmpPath)
    out = "%s/%s.%s.xml" % (tmpPath, genomename, genename)  # Changed from dictionary to tuple
    # Return the parsed variables
    return path, gene, genename, genomename, out


def blastparse(blast_handle, genome, gene, analysisType, cutoff):
    """Parses BLAST results, and populates a dictionary with the results"""
    global plusdict
    records = NCBIXML.parse(blast_handle)   # Open record from memory-mapped file
    dotter()
    incomplete = []
    genomeName = os.path.basename(genome).split('.')[0]
    for record in records:  # This process is just to retrieve HSPs from xml files
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                threadlock.acquire()  # precaution
                # Calculate the percent identity
                percentIdentity = "%.2f" % float(float(hsp.identities) / float(alignment.length) * 100)
                # If the results are greater than the cutoff value, add them to the dictionary
                if hsp.identities >= alignment.length * cutoff:
                    plusdict[genomeName][gene][analysisType] = percentIdentity
                    threadlock.release()
                    # Exit the loop - I added this as else statement was adding
                    # partial matches to the "good" matches found above
                    break
                else:
                    # Make sure that the gene is not already in the dictionary -
                    # may be redundant with the break statement above2
                    if gene not in plusdict[genomeName]:
                        # Puts the HSP in the correct order -  hits to the negative strand will be
                        # reversed compared to what we're looking for
                        if hsp.sbjct_start < hsp.sbjct_end:
                            # Append the start coordinates, end coordinates, and the calculated percent ID
                            incomplete.append((hsp.sbjct_start, hsp.sbjct_end, percentIdentity))
                        else:
                            # Reverse the start and end as required
                            incomplete.append((hsp.sbjct_end, hsp.sbjct_start, percentIdentity))
                    threadlock.release()
    # Once the list is completely populated, find if any partial matches add together for a match that passes the cutoff
    if incomplete:
        # Initialise totalPercentID as the first percent ID value in the list
        totalPercentID = float(sorted(incomplete)[0][2])
        # Initialise adjusted percent ID to 0
        adjPercentID = 0
        # I'm not sure if this is necessary here, but I'm not changing it now
        threadlock.acquire()
        # Initialise currentEntry to the current entry
        currentEntry = [sorted(incomplete)[0][0], sorted(incomplete)[0][1], sorted(incomplete)[0][2]]
        # For each entry in the sorted list of incomplete matches
        # Should look something like: [(1, 915, 45.15), (892, 2048, 49.28)]
        for entry in sorted(incomplete):
            # If entry[0] is less than currentEntry[1], which is less than entry[1]
            # From the example above: 892 is less than 915 is less than 2048
            if entry[0] < currentEntry[1] < entry[1]:
                # The fragment length is the length between entry[0] and entry[1]
                # e.g. 2048 - 892 = 1156
                fragLength = len(range(entry[0], entry[1]))
                # Adjusted fragment length is the length between currentEntry[1] and entry[1]
                # e.g. 2048 - 915 = 1133
                adjFragLength = len(range(currentEntry[1], entry[1]))
                # Adjusted percent ID is the currentEntry percent ID + entry percent ID
                # times the correction factor of adjust fragment length over fragment length
                # e.g. 45.15 + (49.28 * (1133/1156) = 45.15 + (49.28 * 0.9801) = 45.15 + 48.2995 = 93.45
                adjPercentID = float(currentEntry[2]) + (float(entry[2]) * adjFragLength / fragLength)
                # Set current entry to currentEntry[0], entry[1], adjusted percent ID
                # e.g. [1, 2048, 93.45]
                currentEntry = [currentEntry[0], entry[1], "%.2f" % adjPercentID]
                # Update totalPercentID
                totalPercentID = "%.2f" % adjPercentID
            # If entry[0] is greater than currentEntry[1] - simpler calculation,
            # as I don't have to calculate adjusted fragement lengths, I just need to add the percent IDs
            # e.g. [(1, 892, 45.15), (915, 2048, 49.28)] - 915 is greater than 892
            elif entry[0] >= currentEntry[1]:
                # Get the percent ID from entry[2]
                # e.g. 49.28
                adjPercentID = float("%.2f" % float(entry[2]))
                # Add the adjusted percent ID to the currentEntry percentID
                # e.g. 45.15 + 49.28 = 94.43
                totalPercentID = "%.2f" % (float(currentEntry[2]) + adjPercentID)
                # Set the current entry as above
                currentEntry = [currentEntry[0], entry[1], totalPercentID]
        # If the total percent ID calculated above is greater than the cutoff, add the results to the dictionary
        if totalPercentID > cutoff * 100:
            plusdict[genomeName][gene][analysisType] = totalPercentID
        threadlock.release()  # precaution for populate dictionary with GIL


class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)

    def run(self):
        while True:
            global blastpath, plusdict  # global varibles, might be a better way to pipe information
            genome, fasta, blastexist, analysisType, cutoff = self.blastqueue.get()  # retrieve variables from queue
            path, gene, genename, genomename, out = xmlout(fasta, genome)  # retrieve from string splitter
            #Precaution
            threadlock.acquire()
            # Add the appropriate variables to blast path
            blastpath.append((out, path[-1], gene, genename,))  # tuple-list
            try:
                plusdict[genomename][genename] = {analysisType: 0}
            except KeyError:
                plusdict[genomename] = {}
                plusdict[genomename][genename] = {analysisType: 0}
            threadlock.release()
            # Checks to see if this BLAST search has previously been performed
            if not os.path.isfile(out):
                # Print a dot for each gene, genome combination processed
                dotter()
                # Open the output file for writing
                file = open(out, 'w')
                # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db,
                # a mild evalue of 0.1, and XML formatted output
                # Removed perc_identity=percentIdentity from the call, as this allows more flexibility for parsing files
                # with different cutoff values - if I want to re-analyse a search with a lower cutoff, I won't have to
                # re-perform the BLAST search each time
                db = fasta.split(".")[0]
                blastn = NcbiblastnCommandline(query=genome, db=db, evalue=0.1, outfmt=5)
                # Note that there is no output file specified -  the search results are currently stored in stdout
                stdout, stderr = blastn()
                # Search stdout for matches - if the term Hsp appears (the .find function will NOT
                # return -1), a match has been found, and stdout is written to file
                if stdout.find('Hsp') != -1:
                    blast_handle = StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
                    blastparse(blast_handle, genome, genename, analysisType, cutoff)  # parse the data already in memory
                    file.write(stdout)  # write the result
                # Close the file
                file.close()
            # If the BLAST results file already exists is not empty, then parse the results
            elif os.path.getsize(out) != 0:
                # Open the file
                handle = open(out)
                # Read the file into memory
                mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
                # Parse the file in a multithreaded manner
                parsequeue.put((out, genome, genename, mm, analysisType, cutoff))
            # Join all the threads
            parsequeue.join()
            # Error catching?
            if not any(blastpath):
                print out
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes, analysisType, cutoff):
    """Setup and create  threads for blastn and xml path"""
    blastexist = {}
    # Create threads for each gene, genome combination
    for genome in genomes:
        for fasta in fastas:
            # Add the appropriate variables to the threads
            blastqueue.put((genome, fasta, blastexist, analysisType, cutoff))
        blastqueue.join()


class multiparser(threading.Thread): # Had to convert this to a class to integrate threading lock
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)

    def run(self):
        while True:  # General Loop
            global plusdict, genedict  # Import global elements to populate, there may be a better way to do this
            xml, genome, gene, mm, analysisType, cutoff = self.parsequeue.get()  # Retrieve dara from queue
            blastparse(mm, genome, gene, analysisType, cutoff)
            mm.close()
            self.parsequeue.task_done()


def organismChooser(path,targetPath):
    """Allows the user to choose which organism to be used in the analyses"""
    # Initialise a count variable to be used in extracting the desired entry from a list of organisms
    count = 0
    # Check to see if the supplied targetPath has .fa files - if it does, then the default directory structure is probably
    # not being followed, so the target files will be in targetPath
    foldertest = glob("%s/*.fa*" % targetPath)
    if foldertest:
        # Set the required variables as necessary
        # queryGenes are presumably the genes found by foldertest
        queryGenes = foldertest
        # There are likely not going to be qualityGenes included in a custom analysis
        qualityGenes = []
        # Organism name is simply the name of the folder containing the targets
        organismName = targetPath.split("/")[-2]
    else:
        # Get a list of the organisms in the (default) Organism subfolder
        orgList = glob("%sOrganism/*" % path)
        # Iterate through the sorted list
        for folder in sorted(orgList):
            # Ensure that folder is, in actuality, a folder
            if os.path.isdir(folder):
                # Print out the folder names and the count
                print "[%s]: %s" % (count, os.path.split(folder)[1])
                count += 1
        # Get the user input - the number entered corresponds to the list index
        response = input("Please select an organism: ")
        # Get the organism path into a variable
        organism = sorted(orgList)[int(response)]
        organismName = os.path.split(organism)[1]
        # Put the query and quality genes into lists
        queryGenes = glob("%s/query_genes/*.fa" % organism)
        qualityGenes = glob("%s/qualityTest/*.tfa" % organism)
    return queryGenes, qualityGenes, organismName


def blaster(path, cutoff, sequencePath, targetPath):
    """
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    ALL PATHS REQUIRE TRAILING SLASHES!!!
    """
    # Time is used to calculate length of the analyses
    start = time.time()
    # Import global variables
    global count, genedict, blastpath, plusdict
    # Initialise genedict
    genedict = defaultdict(list)
    blastpath = []
    # Run organism chooser to allow the user to choose which databases to use
    # returns the organism name, and lists of
    queryGenes, qualityGenes, organismName = organismChooser(path, targetPath)
    # Get the genome files into a list - note that they must be in the "sequences" subfolder of the path,
    # and the must have a file extension beginning with ".fa"
    strains = glob("%s*.fa*" % sequencePath)
    # Create the threads for the BLAST analysis
    for i in range(len(strains)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        parthreads = multiparser(parsequeue)
        parthreads.setDaemon(True)
        threads.start()
        parthreads.start()
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    # Push targets to threads
    makedbthreads(queryGenes)
    # Quality test genes are optional, so only run the qualityGenes if the folder exists
    if qualityGenes:
        makedbthreads(qualityGenes)
    # There appears to be a race condition going on with the creation of the BLAST databases and the running of BLAST.
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    print "[%s] Now performing and parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    sys.stdout.write('[%s] ' % (time.strftime("%H:%M:%S")))
    # Make blastn threads and retrieve xml file locations
    blastnthreads(queryGenes, strains, "query", cutoff)
    # qualityGenes optional
    if qualityGenes:
        blastnthreads(qualityGenes, strains, "quality", cutoff)
    # Initialise types dictionary
    types = {}
    # Populate types
    types["query"] = queryGenes
    if qualityGenes:
        types["quality"] = qualityGenes
    # print json.dumps(plusdict, sort_keys=True, indent=4, separators=(',', ': '))
    csvheader = ''
    # Loop through the analysis types, and make outputs as required
    for analysisType in types:
        # Initialise variables
        row = ""
        rowcount = 0
        # plusdict contains all results - format: plusdict[genomeName][gene][analysisType] = totalPercentID
        for genomerow in plusdict:
            # The first cell contains the word "Strain"
            csvheader = 'Strain'
            # Append the genome names taken from genomerow
            row += "\n" + genomerow.split('/')[-1].split('.')[0]
            # Increment rowcount
            rowcount += 1
            # For each gene in the appropriate analysis type
            for generow in sorted(types[analysisType]):
                # Extract the gene name from generow
                genename = os.path.basename(generow).split('.')[0]
                # If the gene name is not already in the csvheader variable
                if genename not in csvheader:
                    # Append the gene name in a comma-separated format
                    csvheader += ',' + genename
                    # Format the results - if the percent ID is less than 100
                    identity = plusdict[genomerow][genename][analysisType]
                    if cutoff * 100 < float(identity) < 100:
                    #  Append the percent ID (and a "%") to the row variable
                        row += ',' + str(plusdict[genomerow][genename][analysisType]) + "%"
                    # Otherwise, add a "+" to represent a 100% match
                    elif identity == 0:
                        row += ',N'
                    # If the analysisType is 0 in the dictionary, then there were no matches.
                    # This is shown by an 'N'
                    else:
                        row += ',+'
        # Open the csv report in the appropriate location - add the organism name and the date to keep reports unique
        make_path("%sreports" % path)
        with open("%sreports/%s_%s_results_%s.csv" % (path, organismName, analysisType, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
            # Write the header and the rows
            csvfile.write(csvheader)
            csvfile.write(row)
    # Calculate the elapsed time
    end = time.time() - start
    # Friendly exit statement
    print "\n[%s] Elapsed time for GeneSeeking is %.2f seconds with %.2f seconds per genome" \
          % (time.strftime("%H:%M:%S"), end, end/float(len(strains)))

# Run the blaster function
blaster(path, cutoff, sequencePath, targetPath)
