
import math
import sys
from bunch import Bunch
from fieldtracker import FieldTracker
from orbit.parsers.field_parser import Field_Parser3D


print "Start."


#------------------------------
#Main Bunch init
#------------------------------

b = Bunch()
print "Read Bunch."
runName = "Test_FieldTracker"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.readBunch("parts.dat")
b.getSyncParticle().kinEnergy(energy)
# b.readBunch("Bunches/controlbunch_600.dat", 1)



#------------------------------------------
#Initial Variables
#------------------------------------------
LA11RB  = 0.940222788932
LC11F1  = 0.307193
LC11C12 = 1.8144
LA12RB  = 0.908319544743
XA11M   = 92.55
YA11M   = 23.0
XFOIL1  = 150.4732
YFOIL1  = 46.0
XPFOIL1  = -0.65387068

XMAG = XA11M
YMAG = YA11M

XSTART = XFOIL1
YSTART = YFOIL1
XPSTART = XPFOIL1
YPSTART = 0.0

XREFI = XSTART - XMAG
YREFI = YSTART - YMAG
XREFITOT = 0.0
YREFITOT = 0.0

EULERA = 0.0
EULERB = -0.65387059
EULERG = 0.0
EULERATOT = 0.0
EULERBTOT = 0.0
EULERGTOT = 0.0

ZSTARTTOT = -2.5
ZSTART = LC11F1
ZTOT = 2.5
ZINT = LC11C12 / 2.0
LTOT = ZTOT - ZSTARTTOT
LINT = ZINT - ZSTART
LINTTOT = ZTOT - ZSTART

XPARSEMIN = -21.0
XPARSEMAX =  21.0
YPARSEMIN = -13.0
YPARSEMAX =  13.0

ZPARSEMIN =  100.0 * ZSTART - 1.0
ZPARSEMAX =  100.0 * ZINT + 1.0

SKIPX = 1 
SKIPY = 1
SKIPZ = 1

#=====track bunch through Foil============
 
zsymmetry = 1
printPaths = 1;


mytracker = FieldTracker( 10.190, 10.758, 0.047, 0.056, -0.001, -0.012,
          LINT,  
          ZSTART, ZINT, 0.0001, 2, 1.e-06,
          XREFI, YREFI, EULERA, EULERB, EULERG,b, "testfile.data")

mytracker.setPathVariable(1);

mytracker.trackBunch(b)


b.dumpBunch("final.dat")

print "\nStop."



