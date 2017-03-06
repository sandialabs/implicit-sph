#!/usr/bin/python

import sys, getopt

from itertools import izip_longest

def chunker(iterable, n, fillvalue=None):
   args = [iter(iterable)]*n;
   return izip_longest(*args, fillvalue=fillvalue)

def main(argv):
   mydump  = ''
   mybegin = ''
   myend   = ''
   myout   = ''
   try:
      opts, args = getopt.getopt(argv,"hi:b:e:o:",["input=","begin=","end=","output="])
   except getopt.GetoptError:
      print 'extract.py -i <dump filename> -b <begin of snapshot> -e <end of snapshot> -o <output dump filename>'
      sys.exit(2)

   for opt, arg in opts:
      if opt == '-h':
         print 'extract.py -i <dump filename> -b <begin of snapshot> -e <end of snapshot> -o <output dump filename>'
         sys.exit()
      elif opt in ("-i", "--input"):
         mydump = arg
      elif opt in ("-b", "--begin"):
         mybegin = int(arg)
      elif opt in ("-e", "--end"):
         myend = int(arg)
      elif opt in ("-o", "--output"):
         myout = arg

   timestep = 0;
   natoms = 0;
         
   patt = 'ITEM: TIMESTEP'

   fin  = open(mydump,'r')
   for line in fin:
      if line.startswith(patt):
         timestep = map(int, fin.next().split())[0];
         fin.next();
         natoms = map(int, fin.next().split())[0];         
         break;
   fin.seek(0,0);

   snapsize  = 9 + natoms;
   snapbegin = snapsize*mybegin;
   snapend   = snapsize*(myend+1);

   print "range = ", mybegin, myend
   print "snapshot info = ", snapsize, snapbegin, snapend

   fout = open(myout,'w')   
   for i, line in enumerate(fin):
      if snapbegin <= i < snapend :
         if not i%snapsize:
            print i/snapsize,;
         fout.write(line.rstrip()+'\n');
   fout.close();
   fin.close();
   
if __name__ == "__main__":
   main(sys.argv[1:])
