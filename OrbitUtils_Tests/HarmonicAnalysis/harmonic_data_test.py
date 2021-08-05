#-------------------------------------------------------------------------------
# Class HarmonicData keeps the data y(x) for set of 
# x(i), i=1,Ndata values where each x(i) is a phase between -180 to +180 deg.
# It also has parameters of the fitting representation of this y(x) data:
# y_fit(x) = A0 + A1*cos((PI/180.)*(x + phase_shift_1)) + 
#                 A2*cos((PI/180.)*(2*x + parse_shift_2)) +
#                 ...
# The fast calculation of the sumDiff2 is the main purpose of this class
# This class is used for the harmonics analysis of the scan data for
# the RF cavities when we measure the phase response from BPMs.
#-------------------------------------------------------------------------------

import sys
import math

from orbit_utils import Function
from orbit_utils import HarmonicData

f = Function()

for ind in range(0,360,30):
	x = 1.0*ind
	y = 2.0*math.cos((math.pi/180.)*(x+ 25.))+ 4.0*math.cos((math.pi/180.)*(4*x+ 35.)) + 0.5
	y_err = 0.01*abs(y)
	f.add(x,y,y_err)

order = 4

harmonic_data = HarmonicData(order,f)

harmonic_data.parameter(0,0.5)
harmonic_data.parameter(1,2.0)
harmonic_data.parameter(2,25.0)
harmonic_data.parameter(7,4.0)
harmonic_data.parameter(8,35.0)

print "========================================"
diff2 = harmonic_data.sumDiff2()
print "diff2 =",diff2
print "========================================"

for ind in range(harmonic_data.dataSize()):
	x = harmonic_data.valueX(ind)
	y = harmonic_data.valueY(ind)
	y_err = harmonic_data.valueErr(ind)
	y_fit = harmonic_data.fitValueY(x)
	print "ind= %3d "%ind, " (x,y+-y_err,y_fit) = ( %8.1f , %+8.5f +- %8.5f ,%8.5f)"%(x,y,y_err,y_fit)

print "======================New Data Function ====================="

f.clean()
f.add(0.,1.,0.2)
f.add(1.,3.,0.2)
f.add(2.,4.,0.2)
harmonic_data.setDataFunction(f)

for ind in range(harmonic_data.dataSize()):
	x = harmonic_data.valueX(ind)
	y = harmonic_data.valueY(ind)
	y_err = harmonic_data.valueErr(ind)
	y_fit = harmonic_data.fitValueY(x)
	print "ind= %3d "%ind, " (x,y+-y_err,y_fit) = ( %8.1f , %+8.5f +- %8.5f ,%8.5f)"%(x,y,y_err,y_fit)
	
	
sys.exit(0)

#--------------------------------------------
# This part will check the memory leak 
# in the constructor and resizing the order
# Start the script and open >top command
# in another terminal
#--------------------------------------------

count = 0
while(1<2):
	harmonic_data_new = HarmonicData(harmonic_data)
	count += 1
	harmonic_data_new.order(count%5)
	harmonic_data_new.clean()
	if(count % 10000 == 0):
		print "count = ",count