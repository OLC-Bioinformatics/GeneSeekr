try:
    from setuptools import setup, Command
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.install import install
import sys

import os


class build_card(install):
    description = 'download CARD fasta and modify for ARMI'
    def run(self):
        from bin import makedb
        db = os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'genes.dat')
        if not os.path.isfile(db):
            makedb(db)
        # Attempt to detect whether we were called from setup() or by another
        # command.  If we were called by setup(), our caller will be the
        # 'run_command' method in 'distutils.dist', and *its* caller will be
        # the 'run_commands' method.  If we were called any other way, our
        # immediate caller *might* be 'run_command', but it won't have been
        # called by 'run_commands'.  This is slightly kludgy, but seems to
        # work.
        #
        caller = sys._getframe(2)
        caller_module = caller.f_globals.get('__name__', '')
        caller_name = caller.f_code.co_name
        if caller_module != 'distutils.dist' or caller_name != 'run_commands':
            # We weren't called from the command line or setup(), so we
            # should run in backward-compatibility mode to support bdist_*
            # commands.
            install.run(self)
        else:
            self.do_egg_install()

class updatedb(build_card):
    description = 'update CARD ontology'
    def run(self):
        from bin import updatearo
        updatearo(os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'aro.dat'))
        build_card.run(self)


setup(
    name='pythonGeneSeekr',
    version='0.5dev2',
    packages=['GeneSeekr'],
    package_data={'': ['GeneSeekr/data/aro.dat' ,'GeneSeekr/data/genes.dat', 'GeneSeekr/data/ardb.dat']},
    include_package_data=True,
    url='https://github.com/OLC-Bioinformatics/pythonGeneSeekr',
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='BLAST formatter for full genes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0'],
    cmdclass=dict(install=build_card, card=updatedb),
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
