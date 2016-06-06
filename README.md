# PrimSort
Sort fasta files into groups by primers

##Requirements
Primsort uses :

BioPython  - <https://github.com/biopython/biopython.github.io/>

regex - <https://pypi.python.org/pypi/regex>

##Usage
copy prim_sort2.py, prim_sort.py, fastasplit.py and degreplace.py to same catalogue

**`./prim_sort2.py infasta oligofile no_mismatches_allowed no_processes`**

Will output multiple fasta files named after unique primer names in oligo file.
The fasta files will containe sequences that are a match for primer pairs following
primer name in oligo file.
Also files named infasta.NOMATCH.fa and infasta.AMBIG.fa are generated
containing sequences that do not match any primers or have equally good matches to many.



