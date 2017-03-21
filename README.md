EAModules
==============
A package for modules identification using evolutionary algorithms.

Code structure
==============
- EAModules: core functions of evolutionary algorithms for module identification.
	- functions: core functions of evolutionary algorithms including SA, GA and MA.
	- utils: utilities for supporting core functions like connected component finding and fitness.
	- R: R code of genetic algorithm, modified from COSINE package.
- data: benchmark and real-world data.
- examples: example scripts to call functions, including supplementary files used in the paper.

To construct PPI network, check another package https://github.com/fairmiracle/PPINet.

Reference
==============
[Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. *BMC Genomics* 2017 18(Suppl 2):209](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3495-y)

For any questions, please contact Dong Li at dxl466@cs.bham.ac.uk.
