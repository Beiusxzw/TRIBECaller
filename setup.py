# -*- coding: utf-8 -*-

import re

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.
__version__ = '%s'
"""


def get_version():
    try:
        f = open("TRIBECaller/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)


class install(_install):

    def run(self):
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return

setup(
    name='TRIBECaller',
    version=get_version(),
    author='Ziwei Xue',
    author_email='xueziweisz@gmail.com',
    packages=find_packages(),
    scripts=['bin/CallEditingSites'],
    include_package_data=True,
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        "numpy >= 1.9.0",
        "scipy >= 0.17.0",
        "pysam >= 0.14.0",
        "tqdm >= 4.36.1"
    ],
    entry_points={"console_scripts":["main.py"]}
)