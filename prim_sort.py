#!/usr/bin/python

import sys
import regex
from degreplace import deg_replace
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter

#list to store primer tuples (name,Fprim,Rprim)
primer_list=[]
#populate from oligo file
with open(str(sys.argv[2]), 'r') as f:
 for line in f:
  str_list=line.split('\t')
  primer_list.append((str_list[0],Seq(str_list[1]),Seq(str_list[2].rstrip())))


#list to store search string tuples
#(name,plusFsearch,plusRsearch,minFSearch,minRsearch)
prim_search_strings=[]

#generate and store search string tuples
for primtup in primer_list:
 plusFsearch="(?b)("+deg_replace(str(primtup[1]))+")"+"{e<=2}"
 plusRsearch="(?b)("+deg_replace(str(primtup[2].reverse_complement()))+")"+"{e<=2}"
 minFsearch="(?b)("+deg_replace(str(primtup[1].reverse_complement()))+")"+"{e<=2}"
 minRsearch="(?b)("+deg_replace(str(primtup[2]))+")"+"{e<=2}"
 t=primtup[0],plusFsearch,plusRsearch,minFsearch,minRsearch
 prim_search_strings.append(t)

#create a dictionary  to store match scores
valpairs={}
#create a groupfile
out = open(str(sys.argv[1])+".groups", 'w')

#search sequences for matches
record_iterator = SeqIO.parse(str(sys.argv[1]), "fasta")
#for each sequence
for seq_rec in record_iterator :
#match all posprimers
 for search_tup in prim_search_strings:
  f=regex.search(search_tup[1],str(seq_rec.seq))
  if f is None:
   match_score="NA"
  else:
   r=regex.search(search_tup[2],str(seq_rec.seq))
   if r is None:
    match_score="NA"
   else:
    match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
  valpairs[str(search_tup[0])+"pos"]=match_score
#match all minprimers
 for search_tup in prim_search_strings:
  f=regex.search(search_tup[3],str(seq_rec.seq))
  if f is None:
   match_score="NA"
  else:
   r=regex.search(search_tup[4],str(seq_rec.seq))
   if r is None:
    match_score="NA"
   else:
    match_score=str(sum(f.fuzzy_counts)+sum(r.fuzzy_counts))
  valpairs[str(search_tup[0])+"min"]=match_score
#sort scores
 sortedlist =sorted(valpairs.items(), key=itemgetter(1))
#write line to groupfile
 if sortedlist[0][1] == "NA":
  out.write(seq_rec.id+"\t"+"Nomatch"+"\n")
 else:  
  if sortedlist[0][1]==sortedlist[1][1]:
   out.write(seq_rec.id+"\t"+"Ambig"+"\n")
  else:
   out.write(seq_rec.id+"\t"+sortedlist[0][0]+"\n")

#close group file handle
out.close()








