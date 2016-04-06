# pythonGeneSeekr

### Requirements
The geneSeekr has very few requirements:
 
 1. Python
 2. Sequence (genome) files
 3. Target files
 4. BLAST
 
### Installation
After cloning the git:

```commandline
git clone  https://github.com/OLC-Bioinformatics/GeneSeekr.git
```

Install the python package:

```commandline
python GeneSeekr/setup.py install
```

The following files will be executable from anywhere using commandline

* [GeneSeekr] -- General quick BLASTn formatter 
* [ARMI] -- Antibiotic Resistance Marker Identifier: Use to find markers for any bacterial genome using CARD
* [MLSTSeekr] -- Use with alleles from pubmlst


Upon installation, the `setup.py` script will download the CARD genes to the data folder and will update the CARD onotogies.

Using the `--card_version` and `--card_url` flags you can modify the version of the CARD database to download

Currently `--card_version` defaults to **1.0.4**

*e.g.,*
```commandline
python GeneSeekr/setup.py install --card_version=1.0.4
python GeneSeekr/setup.py install --card_url='some url starting with http://'
```

The `--card_url` flag will trump the `--card_version` flag so it should only used as a last resort

You can bypass the CARD download altogether with which will install the programs into your path
```commandline
python GeneSeekr/setup.py nocard
```

### Examples (Standard Installation)

For [ARMI], simply specify a folder

```commandline
ARMI /path/containing/genomes/in/fasta/format
```

Similarly, for [GeneSeekr] the folder is specfied, along with genes using the `-m` flag (typically inside a single file)

```commandline
GeneSeekr -m ./markers.fasta /path/containing/genomes/in/fasta/format
```

[MLSTSeekr] is a special case of [GeneSeekr] except the markers file is from PubMLST and the names are delimited by an underscore

```commandline
MLSTSeekr  -m ./markers.fasta /path/containing/genomes/in/fasta/format
```


### Commandline Reference
The additional flags can be found using the `--help` flag for the respective program
```commandline
ARMI --help
GeneSeekr --help
MLSTSeekr --help
```

[GeneSeekr]: ./bin/GeneSeekr
[ARMI]: ./bin/ARMI
[MLSTSeekr]: ./bin/MLSTSeekr
