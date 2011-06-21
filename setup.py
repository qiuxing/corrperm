#!/usr/bin/env python

from distutils.core import setup, Extension

corrcperm_module = Extension('_corrcperm',
                             sources=['corrcperm_wrap.cxx', 'corrcperm.cpp'],
                             )

corrcpermtest_module = Extension('corrcpermtest',
                                 sources=['corrcperm.cpp', 'corrcpermtest.cpp'],
                                 )

setup (name = 'corrcperm',
       version = '0.5',
       author      = "Mark Needham, Rui Hu, Sandhya Dwarkadas, Xing Qiu",
       url = "http://www.urmc.rochester.edu/biostat/people/faculty/hu.cfm",
       description = """Detecting Differentially Associated Genes by N-distance and permutation.""",
       ext_modules = [corrcperm_module, corrcpermtest_module],
       py_modules = ['corrcperm'],
       keywords = ["microarray", "statistics", "permutation", "N-distance"],
classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.6",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Bioinformaticians, Biostatisticians",
        "License :: OSI Approved :: GNU General Public License Version 2 (GPL2)",
        "Operating System :: POSIX",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Statistical Analysis :: Bioinformatics",
        ],
)
