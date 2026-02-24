# Phylogenetically Independent Contrasts Logo

## Description

PhyInC Logo (pronounced Fink Logo) for short is a tool to take .fa and .tree files to create
a phylogenetically conditioned sequence logo, using Felsensteins Phylogenetic Independent
Contrast.

## Install from PyPI.org

```
pip3 install phyinc
```
It is always a good idea to use virtual environments, see below.

## Install from a downloaded GitHub repository

For dependencies, start a virtual env as good practice:
```
> python3 -m venv .
> source venv/bin/activate
> pip3 install biopython weblogo matplotlib numpy
```

If you are on mac
```
> brew install ghostscript
```

To install the current package:
```
> pip3 install .
```

## Examples

Run this command:
``` 
> phyinc ./examples/synthetic_data/ex1_t1.tree ./examples/synthetic_data/ex1.fa
```
This should create a PDF named "ex1.fa_seqlogo.pdf" in the examples folder.


## Files Inside

README.md
src/
|-- phyinc/
    |-- config.py  
    |-- phyinc.py  /
    |-- io.py
	|-- version.py
examples/
|-- ex1_t1.tree
|-- ex1.fa


config.py is a configuration file used to save the confirguration from phyinc

ex1_t1.tree is a .tree file provided as an example

ex1.fa is a .fa file provided as an example

phyinc.py is the program

Regular_logo.pdf is the output if you run the program on the example files provided

With_PIC_logo.png is the output if you run the program on the example files provided

# License 

GPLv3


# Authors

* Haolin Guo wrote the basic code as part of his BSc project.
* Kyle Tenn helped make the code into a Python/PyPI package.
* Lars Arvestad oversaw and helped with creating the final Python/PyPI package.

