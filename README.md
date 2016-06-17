# PrimSort
Sort fasta files into groups by primers

##Requirements
Primsort uses :

BioPython  - <https://github.com/biopython/biopython.github.io/>

regex - <https://pypi.python.org/pypi/regex>

##Usage
copy prim_sort2.py, prim_sort.py, fastasplit.py and degreplace.py to same catalogue

**`./prim_sort2.py infasta oligofile mismatches_allowed_nr nr_of_processes`**

Will output multiple fasta files named after unique primer names in oligo file.
The fasta files will containe sequences that are a match for primer pairs following
primer name in oligo file.
Also files named infasta.NOMATCH.fa and infasta.AMBIG.fa are generated
containing sequences that do not match any primers or have equally good matches to many.
Oligo file is a tab separated file with format.
primer_pair_name tab forward_primer tab reverse_primer

##Example
`./prim_sort2.py infile.fa oligos.tsv 2 6`

will sort  infile.fa to several fasta files according to primer_pairs in oligos.tsv
for each primer there can be two mismatches, 6 processes will be started in parallel.
(If you have at least 6 cores/processors the speed will be improved about 6 times)



