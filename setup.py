try:
    from setuptools import setup, Command
    from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.bdist import bdist as _bdist_egg
from distutils.command.build import build as build_module
from GeneSeekr.data import makedb
import os


class bdist_egg(_bdist_egg):
    def run(self):
        self.run_command('Build')
        _bdist_egg.run(self)

class build_card(Command):
    description = 'download CARD fasta and modify for ARMI'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        db = os.path.join(os.path.split(__file__)[0], 'Geneseekr', 'data', 'ARMI-genes.dat')
        if not db:
            makedb(db)
        build_module.build.run(self)

class build(build_module):
    sub_commands = build_module.sub_commands + [('build_card', None)]

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
    cmdclass={'bdist_egg': bdist_egg,
              'build': build,
              'build_card': build_card},
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
