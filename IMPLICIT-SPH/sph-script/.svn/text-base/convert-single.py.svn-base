#!/usr/bin/python

import sys, getopt, re
from dump import dump
from ensight import ensight

def main(argv):
   mydump = ''
   myparaview = ''
   mysingle=''
   try:
      opts, args = getopt.getopt(argv,"hi:o:s",["input=","output=","single="])
   except getopt.GetoptError:
      print 'convert.py -i <dump filename> -o <paraview output name> -s <snapshot number>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convert.py -i <dump filename> -o <paraview output name> -s <snapshot number>'
         sys.exit()
      elif opt in ("-i", "--input"):
         mydump = arg
      elif opt in ("-o", "--output"):
         myparaview = arg
      elif opt in ("-s", "--single"):
         mysingle = arg

   patt = 'ITEM: ATOMS(.*)'
   fd = open(mydump,'r')
   for line in fd:
      m = re.search(patt,line)   
      if m:
         after = m.group(1)
         after = after.strip()
         dumplist = after.split(' ')
         break
   fd.close()

   maplist = list(dumplist);
   for i in range(len(dumplist),0,-1):
       maplist.insert(i-1,i);

   tmplist = list(dumplist[5:]);
   print tmplist;

   paralist = list();
   for i in range(0,len(tmplist),1):
      paralist.append(tmplist[i]);
      paralist.append(tmplist[i]);

   print "maplist:  ", maplist

   d = dump(int(mysingle),mydump);

   print "paraview: ", paralist
   
   d.map(*maplist)
   e = ensight(d);

   print "single snapshot: ", str(mysingle).zfill(4);
   e.single(int(mysingle), myparaview, *paralist);
   
if __name__ == "__main__":
   main(sys.argv[1:])
