#!/usr/bin/python

import sys
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
#from http://biopython.org/wiki/Split_large_file

    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

#input filename and the nr of pieces you want your file split to

def fastasplit(fn,nr):
 #batch size calculation
 record_iter = SeqIO.parse(open(fn),"fasta")
 count=sum(1 for _ in record_iter)
 b_nr = int(count + nr - 1) / nr

 #splitting
 record_iter = SeqIO.parse(open(fn),"fasta")
 for i, batch in enumerate(batch_iterator(record_iter, b_nr)):
  filename = fn+"_%i" % (i + 1)+".temp"
  handle = open(filename, "w")
  count = SeqIO.write(batch, handle, "fasta")
  handle.close()
  print("Wrote %i records to %s" % (count, filename))
