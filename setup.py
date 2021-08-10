from __future__ import print_function
from setuptools import setup, Command

import os


# Functions read() and get_version() were copied from Pip package.
# Purpose is to get version info from current package without it
# being installed (which is usually the case when setup.py is run).
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            # __version__ = "0.9"
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


description = 'LOFAR Solution Tool'
long_description = description
if os.path.exists('README.md'):
    with open('README.md') as f:
        long_description = f.read()


class PyTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import sys, subprocess
        errno = subprocess.call([sys.executable, 'tools/losoto_test.py'])
        raise SystemExit(errno)


setup(
    name='losoto',
    version=get_version('losoto/_version.py'),
    url='http://github.com/revoltek/losoto/',
    project_urls={
        "Documentation": "https://revoltek.github.io/losoto/",
        "Source": "https://github.com/revoltek/losoto/"
    },
    author='Francesco de Gasperin',
    author_email='astro@voo.it',
    license='GPL',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 5 - Production/Stable',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    tests_require=['pytest'],
    install_requires=['numpy>=1.9', 'cython', 'numexpr>=2.0', 'tables>=3.4', 'configparser',
                      'scipy', 'matplotlib', 'python-casacore>=3.0', 'progressbar'],
    scripts=['bin/losoto', 'bin/H5parm_split.py',
             'bin/H5parm2parmdb.py', 'bin/parmdb2H5parm.py', 'bin/killMS2H5parm.py',
             'bin/H5parm_collector.py', 'bin/H5parm_copy.py', 'bin/H5parm_interpolator.py'],
    packages=['losoto', 'losoto.operations', 'losoto.progressbar'],
    test_suite='test',
    cmdclass={'test': PyTest},
)
