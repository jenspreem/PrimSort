#!/usr/bin/python

import os
import sys
from multiprocessing import Process
from prim_sort import prim_sort
from fastasplit import fastasplit

infile=str(sys.argv[1])
oligofile=str(sys.argv[2])
fuzzyval=int(sys.argv[3])
pc=int(sys.argv[4])

fastasplit(infile,pc)


for i in range(1,pc+1):
 if __name__ == '__main__':
  p = Process(target=prim_sort, args=(infile+"_"+str(i)+".temp",'oligos',2))
  p.start()
  p.join()

#remove temp infiles
for i in range(1,pc+1):
 os.remove(infile+"_"+str(i)+".temp")

#concatenate and remove temp outfiles
primnames={'NOMATCH','AMBIG'} # These files are always generated no matter the oligofile

#populate  other primer names from oligo file
with open(str(oligofile), 'r') as f:
 for line in f:
  str_list=line.split('\t')
  primnames.add(str_list[0])

for primname in primnames:
#get names to concatenate
 prim_filenames=[]#keep it ordered
 for i in range(1,pc+1):
  prim_filenames.append(infile+"_"+str(i)+".temp."+primname+".fa")
 with open(infile+"."+primname+".fa","w") as outfile:
  for fname in prim_filenames:
       with open(fname) as inf:
            for line in inf:
               outfile.write(line)
 outfile.close()
#delete temporary outfiles generated for primers
 for fname in prim_filenames:
  os.remove(fname)






