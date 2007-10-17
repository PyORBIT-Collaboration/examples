#!/usr/bin/env python

#---------------------------------------------------------
#This test will read SAD file and parse it
#We have rcs_2.601.sad
#---------------------------------------------------------

import sys
from orbit.parsers.sad_parser import SAD_Parser, SAD_LattElement, SAD_LattLine

if( len(sys.argv) != 2 ):
	print "Usage: >python sad_parser_test.py <name of SAD file>"
	print "Example: >python sad_parser_test.py rsc_2.601.sad"
	sys.exit(1)

sad_file = sys.argv[1]

parser = SAD_Parser()
parser.parse(sad_file)


lines = parser.getSAD_Lines()
elems = parser.getSAD_Elements()
variables = parser.getSAD_Variables()

print "================================================"
print "The whole SAD file includes:"
print "Number of lattice accelerator lines     =",len(lines)
print "Number of lattice accelerator elements  =",len(elems)
print "Number of lattice accelerator variables =",len(variables)
print "================================================"

#get SAD lines dictionary
linesDict = parser.getSAD_LinesDict()

ln = linesDict["INJLINEFL0"]
ln_elems = ln.getElements()
print "The line",ln.getName()," includes N elements=",len(ln_elems)

elem_name = "SBMPD00201"
elem = None
for elem0 in ln_elems:
	if(elem0.getName() == elem_name):
		elem = elem0
		break
if(elem != None):
	print "================================================"
	print "Line:",ln.getName()," (key,val) for element ",elem_name
	for key,val in elem.getParameters().iteritems():
		print "key=",key," \t val=",val
	print "================================================"

print "Done."
sys.exit(0)
