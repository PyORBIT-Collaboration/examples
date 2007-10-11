#!/usr/bin/env python

#---------------------------------------------------------
#This test will read LATTICE and TWISS MAD files, parse them
#and generate the TEAPOT.LAT file. This file is an input for
#original version of the ORBIT for Teapot-like tracking engine.
#---------------------------------------------------------

import sys
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement, MAD_LattLine

if( len(sys.argv) != 5 ):
	print "Usage: >python mad_parser_test.py <name of line> ",
	print "<name of MAD Lattice file> <MAD Twiss File>  <ORBIT Teapot output file> "
	print "Example: >python mad_parser_test.py RING LATTICE TWISS TEAPOT.LAT"
	sys.exit(1)

NAME_LINE = sys.argv[1]
mad_file = sys.argv[2]
twiss_file = sys.argv[3]
teapot_file = sys.argv[4]

parser = MAD_Parser()
parser.parse(mad_file)

lines = parser.getMAD_Lines()
elems = parser.getMAD_Elements()
variables = parser.getMAD_Variables()

print "================================================"
print "The whole MAD file includes:"
print "Number of lattice accelerator lines     =",len(lines)
print "Number of lattice accelerator elements  =",len(elems)
print "Number of lattice accelerator variables =",len(variables)
print "================================================"
print "We will use ",NAME_LINE," lattice line."
print "================================================"

#get MAD lines dictionary
linesDic = parser.getMAD_LinesDic()

ring = linesDic[NAME_LINE]
ring_elems = ring.getElements()
ring_length = 0.
for elem in ring_elems:
	if(elem.getParameters().has_key("L")):
		ring_length = ring_length + elem.getParameter("L")
print "The line",ring.getName()," includes N elements=",len(ring_elems)," length=",ring_length

#--------------------------------------------------
#read TWISS file and calculate number of elements
#--------------------------------------------------
twiss_file = open(twiss_file,"r")
i = 0
for line in twiss_file:
	i = i + 1
#3 lines in the begginning
#5 lines - initial point twiss
#1 line ststistics at the end
#2 lines H and V tune, twiss etc.
i = i - 2 - 5 - 1 - 2
if(i%5 != 0):
	print "Twiss file has unusual structure! Take a look at Twiss!"
	sys.exit(1)
i = i / 5
twiss_file.close()
print "Twiss file has ", i," elements."

#-------------------------------------------------
#Let us write down TEAPOT file
#-------------------------------------------------
teapot_file = open(teapot_file,"w")
for elem in ring_elems:
	res = elem.getType()+" "
	for (key,val) in elem.getParameters().items():
		res = res + key + " " + str(val) + " "
	res = res + "end \n"
	teapot_file.write(res.lower())


teapot_file.write("final \n")
teapot_file.close()

print "Done."
sys.exit(0)
