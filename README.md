# pythonGeneSeekr

The geneSeekr has very few requirements:
 
 1. Python
 2. BioPython
 3. Sequence (genome) files
 4. Target files
 5. ?
 
 There are three flags that must be provided when running the program from a system other than my own. These flag override hard-coded paths.

After cloning the repository running `python setup.py install` will install 

* GeneSeekr -- General quick BLASTn formatter 
* ARMI -- Antibiotic Resistance Marker Identifier: Use to find markers for any bacterial genome using CARD
* MLSTSeekr -- Use with alleles from pubmlst

These programs will become a part of the `$PATH`

Upon installation, the `setup.py` script will download the CARD genes to the data folder. Running `setup.py install card` will update the CARD onotogies into the `data/aro.dat` pickle file

```
usage: GeneSeekr [-h] [--version] -i INPUT [-o OUTPUT] [-c CUTOFF]
                 [-t THREADS] -m MARKER

General quick BLASTn formatter

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        Specify input fasta folder
  -o OUTPUT, --output OUTPUT
                        Specify output folder for csv and json
  -c CUTOFF, --cutoff CUTOFF
                        Threshold for maximum unique bacteria for a single
                        allele
  -t THREADS, --threads THREADS
                        Specify number of threads
  -m MARKER, --marker MARKER
                        Specify gene file in FASTA format
```

```
usage: ARMI [-h] [--version] -i INPUT [-o OUTPUT] [-c CUTOFF] [-t THREADS]
            [-m MARKER] [-a ANTI] [--tolc]

Antibiotic Resistance Marker Identifier: Use to find markers for any bacterial
genome

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        Specify input fasta folder
  -o OUTPUT, --output OUTPUT
                        Specify output folder for csv and json
  -c CUTOFF, --cutoff CUTOFF
                        Threshold for maximum unique bacteria for a single
                        allele
  -t THREADS, --threads THREADS
                        Specify number of threads
  -m MARKER, --marker MARKER
                        Specify gene file in FASTA format
  -a ANTI, --anti ANTI  JSON file location
  --tolc                Include TolC-related efflux pumps
```

```
usage: MLSTSeekr [-h] [--version] -i INPUT [-o OUTPUT] [-c CUTOFF]
                 [-t THREADS] -m MARKER

Multilocus Seqeunce Typing Assay with BLAST: Use to find markers for any
bacterial genome

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        Specify input fasta folder
  -o OUTPUT, --output OUTPUT
                        Specify output folder for csv and json
  -c CUTOFF, --cutoff CUTOFF
                        Threshold for maximum unique bacteria for a single
                        allele
  -t THREADS, --threads THREADS
                        Specify number of threads
  -m MARKER, --marker MARKER
                        Specify gene file in FASTA format
```
