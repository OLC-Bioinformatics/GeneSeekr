try:
    from setuptools import setup, Command
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.install import install
from GeneSeekr.data import makedb, updatearo
import os


class build_card(install):
    description = 'download CARD fasta and modify for ARMI'
    def run(self):
        db = os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'genes.dat')
        if not os.path.isfile(db):
            makedb(db)
        install.run(self)

class updatedb(build_card):
    description = 'update CARD ontology'
    def run(self):
        updatearo(os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'aro.dat'))
        build_card.run(self)


setup(
    name='pythonGeneSeekr',
    version='0.5dev2',
    packages=['GeneSeekr'],
    package_data={'': ['GeneSeekr/data/aro.dat','GeneSeekr/data/genes.dat', 'GeneSeekr/data/ardb.dat']},
    include_package_data=True,
    url='https://github.com/OLC-Bioinformatics/pythonGeneSeekr',
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='BLAST formatter for full genes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0'],
    cmdclass={'install': build_card,
              'card': updatedb},
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
