# PrimSort
Sort fasta files into groups by primers/barcodes/oligos etc.

##Requirements
Primsort uses :

BioPython  - <https://github.com/biopython/biopython.github.io/>

regex - <https://pypi.python.org/pypi/regex>

##Usage
copy prim_sort.py and degreplace.py to same catalogue

**`./prim_sort.py fasta_to_sort`**

Will output  **fasta_to\_sort.groups** file containing fasta 
sequence headers (without >) and group names separated by tab


