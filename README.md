# Phylogenetically Independent Contrasts Logo

## PhyInC Logo (pronounced Fink Logo) for short is a tool to take .fa and .tree files to create sequence Logo and save as a .PNG file.

For dependencies
>pip install biopython weblogo matplotlib numpy

If you are on mac
>brew install ghostscript

Run Example
>python3 PIC.py ex1_t1.treefile ex1.fa

## Files Inside

README.md    
Code/  
|-- config.py  
|-- ex1_t1.tree  
|-- ex1.fa  
|-- phyinc.py  
|-- Regular_logo.pdf  
`-- With_PIC_logo  

config.py is a configuration file used to save the confirguration from phyinc

ex1_t1.tree is a .tree file provided as an example

ex1.fa is a .fa file provided as an example

phyinc.py is the program

Regular_logo.pdf is the output if you run the program on the example files provided

With_PIC_logo.png is the output if you run the program on the example files provided

# License 

GPLv3