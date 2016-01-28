try:
    from setuptools import setup, Command
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.install import install
from GeneSeekr.data import makedb
import os


class build_card(install):
    description = 'download CARD fasta and modify for ARMI'
    def run(self):
        db = os.path.abspath(os.path.join(os.path.split(__file__)[0], 'Geneseekr', 'data', 'ARMI-genes.dat'))
        if not os.path.isfile(db):
            makedb(db)
        install.run(self)


setup(
    name='pythonGeneSeekr',
    version='0.5dev2',
    packages=['GeneSeekr'],
    package_data={'': ['GeneSeekr/data/*.dat']},
    include_package_data=True,
    url='https://github.com/OLC-Bioinformatics/pythonGeneSeekr',
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='BLAST formatter for full genes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0'],
    cmdclass={'install': build_card},
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
