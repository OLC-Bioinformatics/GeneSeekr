from distutils.core import setup
from distutils.command import build as build_module
from GeneSeekr.data import makedb
import os


class Build(build_module.build):
    def run(self):
        db = os.path.join(os.path.split(__file__)[0], 'Geneseekr', 'data', 'ARMI-genes.dat')
        if not db:
            makedb(db)
        build_module.build.run(self)

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
    cmdclass={'build': Build},
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
