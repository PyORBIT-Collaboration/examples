#! /usr/bin/env python

"""
This is a test script to check the functionality of the 
linac acc lattice.
"""

import sys

from linac_parser import SimplifiedLinacParser
from LinacAccLattice import LinacLatticeFactory, LinacAccLattice

parser = SimplifiedLinacParser("sns_linac.xml")
linacTree = parser.getLinacStructTree()
print "======================================="
print "Total length=",linacTree.getLength()
print "======================================="
sequences = linacTree.getSeqs()
for seq in sequences:
	print "seq=",seq.getName()," L=",seq.getLength()
	
lattFactory = 	LinacLatticeFactory(linacTree)
accLattice = lattFactory.getLinacAccLattice(["MEBT","DTL1","DTL2"])
	
sys.exit(1)

