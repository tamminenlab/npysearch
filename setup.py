import sys

from pybind11 import get_cmake_dir
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension("_npysearch",
        ["src/_npysearch/Search.cpp", "src/_npysearch/TextReader.cpp"],
        define_macros = [('VERSION_INFO', __version__)],
        include_dirs  = ["include"]
        ),
]

setup(
    name             = "npysearch",
    version          = __version__,
    author           = "Aditya Jeevannavar",
    author_email     = "adjeev@utu.fi",
    url              = "https://github.com/jeevannavar/npysearch",
    description      = "A test project to implement python bindings for nsearch",
    long_description = "",
    packages         = ["npysearch"],
    ext_modules      = ext_modules,
    extras_require   = {"test": "pytest"},
    cmdclass         = {"build_ext": build_ext},
    zip_safe         = False,
    python_requires  = ">=3.6",
)