import sys
from bunch import Bunch

#-----------------------------------------------------
#Test of particle attributes for bunch
#-----------------------------------------------------

print "Start."

b = Bunch()

nParts = 5
for i in xrange(nParts):
	b.addParticle(0.1+i,0.2+i,0.3+i,0.4+i,0.5+i,0.6+i)
b.compress()

# add particles attributes with the optional parameters dictionary 
d = {"size":5}
b.addPartAttr("Amplitudes",d)

b.dumpBunch("bunch_dump_test.dat")

b_new = Bunch()
b_new.readBunch("bunch_dump_test.dat")
b_new.addPartAttr("macrosize")
b_new.dumpBunch()

print "part. attr. dicts=",b_new.getPartAttrDicts()

nms = b_new.readPartAttrNames("bunch_dump_test.dat")
print "debug names of attr. =",nms

res_dict = b_new.readPartAttrDicts("bunch_dump_test.dat")
print "debug dict=",res_dict

print "Stop."

sys.exit(1)

count = 0
while(1 < 2):
	count = count + 1
	b_new = Bunch()
	b_new.readBunch("bunch_dump_test.dat")
	if(count % 100 == 0): print "i=",count
