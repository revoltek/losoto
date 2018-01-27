from __future__ import print_function
from setuptools import setup, Command

import os

import losoto._version

description = 'LOFAR Solution Tool'
long_description = description
if os.path.exists('README.md'):
    with open('README.md') as f:
        long_description=f.read()


class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'tools/losoto_test.py'])
        raise SystemExit(errno)


setup(
    name='losoto',
    version=losoto._version.__version__,
    url='http://github.com/revoltek/losoto/',
    author='Francesco de Gasperin',
    author_email='astro@voo.it',
    description=description,
    long_description=long_description,
    platforms='any',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: Stable',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    tests_require=['pytest'],
    install_requires=['numpy>=1.9','cython','numexpr>=2.0','tables>=3.0'],
    scripts = ['bin/losoto', 'bin/H5parm_benchmark.py',
               'bin/H5parm_exporter.py', 'bin/H5parm_importer.py',
               'bin/H5parm_collector.py','bin/H5parm_copy.py'],
    packages=['losoto','losoto.operations','losoto.progressbar'],
    test_suite='test',
    cmdclass = {'test': PyTest},
    )
