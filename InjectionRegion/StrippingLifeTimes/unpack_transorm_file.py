#!/usr/bin/env python

#---------------------------------------------------------
#This will beautify the transition file for Timofei Gorlov
# personally.
#---------------------------------------------------------

import sys

if( len(sys.argv) != 3 ):
	print "Usage: >python ",sys.argv[0]," transition.txt  beauty_transition.txt"
	sys.exit(1)
	
fl_in = open(sys.argv[1],"r")
fl_out = open(sys.argv[2],"w")

count = 0
for line in fl_in:
	line = line.expandtabs(1)
	if(count%1000 == 0): print "line #=",count
	trans_part = line[:22]
	data_arr = line[23:].split()
	st = trans_part.ljust(25)
	for ind in range(len(data_arr)):
		dt = data_arr[ind].strip()
		if(ind%3 == 0): st = st + "        "
		if(dt == "0.0" ):
			st = st + "          0.000000000000000000000000000000  "
		else:
			np = dt.find(".")
			n = 11 - np
			for ii in range(n): st = st + " "
			st = st + dt + "  "
	fl_out.write(st+"\n")
	count = count + 1

fl_out.close()
fl_in.close()


	
