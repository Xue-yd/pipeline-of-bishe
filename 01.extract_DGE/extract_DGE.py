#!/usr/bin/python

import re
fh1 = open("human_gene.txt")
fh2 = open("GSE142585_MergedHumanTestis4_DGE.txt")
fh3 = open("human_orthologs_DGE.txt","w")

firstline = fh2.readline()
fh3.write(firstline)
fh2.close()
for line1 in fh1:
      name = "^"+line1.strip()+"\s"
      fh2 = open("GSE142585_MergedHumanTestis4_DGE.txt")
      firstline = fh2.readline()
      for line2 in fh2:
         rgx = re.compile(name)
         result = rgx.search(line2)
         if(result != None):
            fh3.write(line2)
      fh2.close()

fh1.close()
fh3.close()
   
