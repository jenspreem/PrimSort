#!/usr/bin/python

import sys
import regex
from degreplace import deg_replace
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter

#bakterite forward
bac_prim=Seq("GTGYCAGCMGCCGCGGTAA")
#bakterite reverse 
bac_rev=Seq("AGTCAGCCAGCCGYCAATTYMTTTRAGTTT")

#arhede forward
arc_prim=Seq("CAGYCGCCRCGGTAA")
#arhede reverse
arc_rev=Seq("GCYCCCCCGCCWATTC")

#generate search strings to search + strand .. deg_replace replaces degenrate nucleotides
#with appropriate nucleotide groups
bacFsearch="(?b)("+deg_replace(str(bac_prim))+")"+"{e<=2}"
bacRsearch="(?b)("+deg_replace(str(bac_rev.reverse_complement()))+")"+"{e<=2}"

arcFsearch="(?b)("+deg_replace(str(arc_prim))+")"+"{e<=2}"
arcRsearch="(?b)("+deg_replace(str(arc_rev.reverse_complement()))+")"+"{e<=2}"

#generate search strings for - strands
bacFminsearch="(?b)("+deg_replace(str(bac_prim.reverse_complement()))+")"+"{e<=2}"
bacRminsearch="(?b)("+deg_replace(str(bac_rev))+")"+"{e<=2}"

arcFminsearch="(?b)("+deg_replace(str(arc_prim.reverse_complement()))+")"+"{e<=2}"
arcRminsearch="(?b)("+deg_replace(str(arc_rev))+")"+"{e<=2}"

#create a dictionary  to store match scores
valpairs={}
#create a groupfile
out = open(str(sys.argv[1])+".groups", 'w')

#search sequences for matches
record_iterator = SeqIO.parse(str(sys.argv[1]), "fasta")
for seq_rec in record_iterator :

#bac + 
 f=regex.search(bacFsearch,str(seq_rec.seq))
 if f is None:
  match_score="NA"
 else:
  r=regex.search(bacRsearch,str(seq_rec.seq))
  if r is None:
   match_score="NA"
  else:
   match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
 valpairs["Bpos"]=match_score

#bac -
 f=regex.search(bacFminsearch,str(seq_rec.seq))
 if f is None:
  match_score="NA"
 else:
  r=regex.search(bacRminsearch,str(seq_rec.seq))
  if r is None:
   match_score="NA"
  else:
   match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
 valpairs["Bmin"]=match_score

#arc+

 f=regex.search(arcFsearch,str(seq_rec.seq))
 if f is None:
  match_score="NA"
 else:
  r=regex.search(arcRsearch,str(seq_rec.seq))
  if r is None:
   match_score="NA"
  else:
   match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
 valpairs["Apos"]=match_score

#arc-
 f=regex.search(arcFminsearch,str(seq_rec.seq))
 if f is None:
  match_score="NA"
 else:
  r=regex.search(arcRminsearch,str(seq_rec.seq))
  if r is None:
   match_score="NA"
  else:
   match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
 valpairs["Amin"]=match_score
#sort scores
 sortedlist =sorted(valpairs.items(), key=itemgetter(1))
#write line to groupfile
 if sortedlist[0][1] != "NA":
  out.write(seq_rec.id+"\t"+sortedlist[0][0]+"\n")
 

#close group file handle
out.close()





