# -*- coding: utf-8 -*-
import os
import codecs
from distutils.core import setup
from distutils.cmd import Command

def readfile(filename):
    with codecs.open(filename,  encoding='utf-8') as f:
        return f.read()

setup(
    name='hecke_hurwitz',
    version='0.0',
    description='A Python module that depends on Sage to discover Hecke-Hurwitz GL(2,R)-orbit closure',
    long_description=readfile("README.rst"),
    author='Vincent Delecroix',
    author_email='vincent.delecroix@u-bordeaux.fr',
    # Remember to update these if the directory structure changes.
    packages=['hecke_hurwitz/'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ]
)
