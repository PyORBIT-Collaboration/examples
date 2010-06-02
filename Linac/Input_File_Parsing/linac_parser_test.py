#! /usr/bin/env python

"""
This is a test script to check the functionality of the 
linac parser.
"""

import sys

from orbit.sns_linac  import SimplifiedLinacParser

parser = SimplifiedLinacParser("sns_linac.xml")
linacTree = parser.getLinacStructureTree()
print "======================================="
print "Total length=",linacTree.getLength()
print "======================================="
sequences = linacTree.getSeqs()
for seq in sequences:
	print "seq=",seq.getName()," L=",seq.getLength()
	
print "======================================="	
nodes = sequences[0].getNodes()
for node in nodes:
	print "node=",node.getName()," type=",node.getType()," position=",node.getParam("pos")," L=",node.getLength()



sys.exit(1)

