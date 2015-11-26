from distutils.core import setup

setup(
    name='pythonGeneSeekr',
    version='0.5dev',
    packages=['GeneSeekr'],
    url='https://github.com/OLC-Bioinformatics/pythonGeneSeekr',
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='BLAST formatter for full genes',
    long_description=open('README.md').read(),
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0'],
    scripts=['bin/ARMI',
             'bin/GeneSeekr']
)
