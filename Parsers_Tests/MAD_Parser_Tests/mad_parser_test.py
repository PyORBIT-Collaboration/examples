#!/usr/bin/env python

#---------------------------------------------------------
#This test will read MAD file and parse it
#We have SNSring.mad (it will call to SNSring.lat)
#---------------------------------------------------------

import sys
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement, MAD_LattLine

if( len(sys.argv) != 2 ):
	print "Usage: >python mad_parser_test.py <name of MAD file>"
	print "Example: >python mad_parser_test.py SNSring.mad"
	sys.exit(1)

mad_file = sys.argv[1]

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

#get MAD lines dictionary
linesDic = parser.getMAD_LinesDic()

ring = linesDic["RNG"]
ring_elems = ring.getElements()
print "The line",ring.getName()," includes N elements=",len(ring_elems)

print "Done."
sys.exit(0)
