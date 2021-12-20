# npysearch

npysearch implements an efficient BLAST-like sequence comparison algorithm, written in C++11 and using native Python datatypes and bindings. npysearch is light-weight, fast, and dependency-free. The code base of npysearch is adapted from [nsearch](https://github.com/stevschmid/nsearch).


## Installation

### from pypi

```
pip install npysearch
```

### from conda-forge

Installing from the `conda-forge` channel can be achieved by adding `conda-forge` to your channes with:
```
conda config --add channels conda-forge
conda config --set channel_priority strict
```
You can skip the above step if `conda-forge` channel has been added already.

Once the conda-forge channel has been enabled, `npysearch` can be installed with:

```
conda install npysearch
```

### from github

```
# Clone repository from github
git clone https://github.com/jeevannavar/npysearch.git

# Install package using pip
pip install ./npysearch
```


## Examples

```Python
# Import npysearch package
import npysearch as npy

# Read query file into a dictionary
query = npy.read_fasta("npysearch/inst/extdata/query.fasta")

# Read database file into a dictionary
database = npy.read_fasta("npysearch/inst/extdata/db.fasta")

# BLAST the query against the database
results_dna = blast(query, database)

# BLAST protein sequence file against itself using filenames as blast function arguments

results_prot = blast(query = "npysearch/inst/extdata/prot.fasta",
                     database = "npysearch/inst/extdata/prot.fasta",
                     alphabet = "protein")
```

## Caveats

* The `blast` function automatically detects whether the query and database arguments were passed as string paths to fasta files or as dictionaries of sequences. Both of them need not be input as the same type.
* Use `help(npy)` (assuming you've imported npysearch as npy) to get a list of the functions included and their docstrings. For docstrings of specific functions, for example blast, use `help(npy.blast)`