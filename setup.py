import sys
import os

from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

__version__ = "1.3.0"

this_dir = os.path.dirname(os.path.abspath(__file__))

ext_modules = [
    Pybind11Extension("_npysearch",
        ["src/_npysearch/Search.cpp", "src/_npysearch/TextReader.cpp"],
        define_macros = [('VERSION_INFO', __version__)],
        include_dirs  = [os.path.join(this_dir, "include"),
                         os.path.join(this_dir, "src/_npysearch")]
        ),
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 5 - Production/Stable',

    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',

    # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: BSD License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.9'
]

setup(
    name             = "npysearch",
    version          = __version__,
    author           = "Aditya Jeevannavar",
    author_email     = "adjeev@utu.fi",
    url              = "https://github.com/jeevannavar/npysearch",
    description      = "Python bindings for nsearch, an efficient BLAST-like sequence comparison algorithm written in C++",
    long_description = long_description,
    long_description_content_type="text/markdown",
    packages         = ["npysearch"],
    ext_modules      = ext_modules,
    extras_require   = {"test": "pytest"},
    cmdclass         = {"build_ext": build_ext},
    zip_safe         = False,
    python_requires  = ">=3.6",
    classifiers      = classifiers,
    keywords         = "BLAST nsearch bioinformatics"
)
