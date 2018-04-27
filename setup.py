import os
import sys

from distutils.core import setup
from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

install_requires = ["pandas"]
# try:
#     os.getenv("TRAVIS")
#     install_requires.append("coveralls")
# except:
#     pass

# if sys.version_info[0] == 2:
#     install_requires.append("functools32")

macros = [("CYTHON_TRACE", "1")]
# macros = Nonec

if macros:
    from Cython.Compiler.Options import get_directive_defaults
    directive_defaults = get_directive_defaults()
    directive_defaults['linetrace'] = True
    directive_defaults['binding'] = True


e1 = Extension("src.intervaltree", ["src/intervaltree.pyx"], define_macros = macros)

extensions = [e1]

setup(
    name="pyranges",
    packages=find_packages(),
    package_data={'pyranges': ['example_data/*.bed']},
    ext_modules=cythonize(extensions),
    include_dirs=["."],
    version="0.0.1",
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
