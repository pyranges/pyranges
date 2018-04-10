import os
import sys
from setuptools import setup, find_packages
# from Cython.Build import cythonize

from pyranges.version import __version__
install_requires = ["pandas"]
# try:
#     os.getenv("TRAVIS")
#     install_requires.append("coveralls")
# except:
#     pass

# if sys.version_info[0] == 2:
#     install_requires.append("functools32")

setup(
    name="pyranges",
    packages=find_packages(),

    # scripts=["bin/featurefetch"],
    version=__version__,
    description="GRanges for Python.",
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
    ("Pythonic Genomic Ranges."))
