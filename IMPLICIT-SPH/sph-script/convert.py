#!/usr/bin/python

import sys, getopt, re
from dump import dump
from ensight import ensight

def main(argv):
   mydump = ''
   myparaview = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["input=","output="])
   except getopt.GetoptError:
      print 'convert.py -i <dump filename> -o <paraview output name>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'convert.py -i <dump filename> -o <paraview output name>'
         sys.exit()
      elif opt in ("-i", "--input"):
         mydump = arg
      elif opt in ("-o", "--output"):
         myparaview = arg

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
   print "paraview: ", paralist

   d = dump(mydump);
   
   d.map(*maplist)
   e = ensight(d);
   e.one(myparaview, *paralist);
   
if __name__ == "__main__":
   main(sys.argv[1:])
