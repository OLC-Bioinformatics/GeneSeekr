try:
    from setuptools import setup, Command
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup, Command
    from distutils.command.install import install
import sys
import os


class UpdateDB(install):
    description = 'update CARD ontology'
    user_options = install.user_options + [
        ('card-version=', None, 'Specify the version of CARD to download (default 1.0.3)'),
        ('card-url=', None, 'URL of CARD to download, this should only be used if card_version doesn\'t work')
    ]

    def run(self, obj=None):
        from GeneSeekr.data import Build
        card = Build(self.card_url)
        card.updatearo(os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'aro.dat'))
        card.makedb(os.path.join(os.path.split(__file__)[0], 'GeneSeekr', 'data', 'genes.dat'))
        del card
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

    def initialize_options(self):
        install.initialize_options(self)
        self.card_version = None
        self.card_url = None

    def finalize_options(self):
        install.finalize_options(self)
        self.card_version = self.card_version if self.card_version else "1.0.4"
        if self.card_version and not self.card_url:
            assert self.card_version != '1.0.0', 'Version 1.0.0 not supported'
            assert "." in self.card_version, 'Invalid Card Version {0!r:s}'
            assert all(map(str.isdigit, self.card_version.split("."))), 'Invalid Card Version {0!r:s}'
            self.card_url = 'https://card.mcmaster.ca/download/0/broadsteet-v{0:s}.tar.gz'.format(self.card_version)
        elif self.card_url:
            assert self.card_url.startswith("http"), "URL must start with http "
        import urllib2
        request = urllib2.Request(self.card_url)
        request.get_method = lambda: 'HEAD'
        try:
            response = urllib2.urlopen(request)
            assert int(response.getcode()) < 400, \
                'HTTP Error {0:s}: ({1:s})'.format(response.response.getcode(), self.card_url)
            assert response.info().type != 'text/html', 'URL {0:s} points to webpage not file'.format(self.card_url)
        except urllib2.HTTPError as e:
            print 'HTTP Error {0:s}: {1:s} ({2:s})'.format(e.code, e.msg, e.url)
            raise
        except:
            print 'Malformed CARD URL Error'
            print u'If you do did not specify the URL:\'{0:s}\' or version:\'{1:s}\' then either CARD is down or ' \
                  u'the link has changed, you may be able to remedy this by specifying a url or version with ' \
                  u'--card-version=\'version\' or  --card-url=\'url\' after \'{2:s}\''\
                .format(self.card_url, self.card_version, ' '.join(sys.argv))
            raise

setup(
    name='GeneSeekr',
    version='0.6.dev1',
    packages=['GeneSeekr'],
    package_data={'': ['GeneSeekr/data/*.dat']},
    include_package_data=True,
    url='https://github.com/OLC-Bioinformatics/GeneSeekr',
    license='MIT',
    author='mike knowles',
    author_email='mikewknowles@gmail.com',
    description='BLAST formatter for full genes',
    long_description=open('README.md').read(),
    setup_requires=['setuptools >= 18.0.1'],
    install_requires=['biopython >= 1.65',
                      'argparse >= 1.4.0',
                      'pysam == 0.8.4',
                      'pysamstats'],
    cmdclass=dict(install=UpdateDB, nocard=install),
    scripts=['bin/ARMI',
             'bin/MLSTSeekr',
             'bin/GeneSeekr',
             'bin/resfinder']
)
