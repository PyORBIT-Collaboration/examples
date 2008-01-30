import sys
from bunch import Bunch

#-----------------------------------------------------
#Dump and read bunch to and from file
#-----------------------------------------------------

print "Start."

b = Bunch()

nParts = 5
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
b.compress()

# add particles attribute
b.addPartAttr("macrosize")

#add bunch attributes
b.bunchAttrDouble("aaa",222)
b.bunchAttrDouble("bbb",333)
b.bunchAttrInt("aaa",111)
b.bunchAttrInt("bbb",444)

b.dumpBunch("bunch_dump_test.dat")

b_new = Bunch()
b_new.readBunch("bunch_dump_test.dat")

b_new.dumpBunch()

print "Stop."
