import os
import sys

from distutils.core import setup
from setuptools import find_packages

install_requires = ["cython", "pandas", "ncls", "clustertree", "tabulate", "pysam", "sorted_nearest"]



setup(
    name="pyranges",
    packages=find_packages(),
    package_data={'pyranges': ['example_data/*.bed', 'example_data/*.gtf',
                               'example_data/*.bam', 'example_data/*.bam.bai']},
    include_dirs=["."],
    version="0.0.4",
    description="PyRanges for Python.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/pyranges",
    keywords=["Bioinformatics"],
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment", "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    long_description=
    ("Performant Pythonic Genomic Ranges."))
