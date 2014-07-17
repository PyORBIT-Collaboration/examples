#!/usr/bin/env python

#--------------------------------------------------------
# The classes will read data file with lattice information
#--------------------------------------------------------

import math
import sys
import os

fl_in = open("trace_win_scl_structure.dat","r")
lns = fl_in.readlines()
fl_in.close()

count = 0
quads = []
for ind in range(len(lns)):
	ln = lns[ind]
	ln = ln.strip()
	if(len(ln) <= 0 or ln[0] == ";"): continue
	if(ln.find("QUAD") >= 0):
		count += 1
		if(lns[ind-1].find("MATCH_FAM_GRAD") >= 0):
			print "quad count=",count," ln=",ln
			
print "total quads N=",count
