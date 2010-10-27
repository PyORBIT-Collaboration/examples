#-----------------------------------------------------
#Creates Space Charge Calculator based on Rick Baartman 
#suggestion and tracks the test bunch through the calculator. 
#-----------------------------------------------------


import sys
import math
import random

import orbit_mpi

from bunch import Bunch
from spacecharge import SpaceChargeCalc2p5Drb

print "Start."

sizeX = 64
sizeY = 64
sizeZ = 5
long_avg_n = 3
calc2p5d = SpaceChargeCalc2p5Drb(sizeX,sizeY,sizeZ)
calc2p5d.setLongAveragingPointsN(long_avg_n)

macroSize = 1.0e+13
energy = 1.0
b = Bunch()
b.macroSize(macroSize)
b.getSyncParticle().kinEnergy(energy)
gamma = b.getSyncParticle().gamma()
beta = b.getSyncParticle().beta()
print "gamma=",b.getSyncParticle().gamma()
print "beta=",b.getSyncParticle().beta()
#b.addPartAttr("macrosize")

bunch_radius = 0.005
bunch_length = 200.0

nParts = 1000000

for ip in range(nParts):
	r = bunch_radius*math.sqrt(random.random())
	phi = 2*math.pi*random.random()
	x = r*math.sin(phi)
	y = r*math.cos(phi)
	z = bunch_length*0.5*(1.0 - 2*random.random())
	"""
	z = 0.5*bunch_length*math.sqrt(random.random())
	if(random.random() > 0.5):
		z = z - bunch_length/2
	else:
		z = bunch_length/2 - z
	"""
	b.addParticle(x,0.,y,0.,z,0.)

print "bunchSize = ",b.getSize()
print "macroSize=",b.macroSize()
print "mass=",b.mass()

pipe_radius = 0.010
slice_length = 0.1

#b.dumpBunch("pyorbit_bunch_test_in.dat")
print "Start Poisson Solver."
calc2p5d.trackBunch(b,slice_length,pipe_radius)
print "Stop Poisson Solver."
#b.dumpBunch("pyorbit_bunch_test_out.dat")

rhoGrid = calc2p5d.getRhoGrid()
phiGrid = calc2p5d.getPhiGrid()
longGrid = calc2p5d.getLongGrid()
longDerivGrid = calc2p5d.getLongDerivativeGrid()

x = bunch_radius/2.0
y = 0.
r = math.sqrt(x*x+y*y)

rho = rhoGrid.getValue(x,y)
rho_theory = b.macroSize()*4.0/(math.pi*rhoGrid.getSizeX()*rhoGrid.getSizeY())
print "r=",r," rho  = %12.5g "%rho,"  rho_theory = %12.5g "%rho_theory

#--------------------------------------------------------------------------
# remember that the potential is an abstract potential: lambda*ln(r)
# The potential for a charged string in CGS is 2*lambda*ln(r)
#--------------------------------------------------------------------------
phi = 2*(phiGrid.getValue(x,y) - phiGrid.getValue(0.,0.))
phi_theory = b.macroSize()*r**2/(bunch_radius**2)
print "r=",r," phi  = %12.5g "%phi,"  phi_theory = %12.5g "%phi_theory

(ex,ey) = phiGrid.calcGradient(x,y)
grad = 2*math.sqrt(ex*ex+ey*ey)
grad_theory = b.macroSize()*2*r/(bunch_radius**2)
print "r=",r," grad = %12.5g "%grad," grad_theory = %12.5g "%grad_theory

# theoretical coeff delta(r_prime/r)
slope_theory = (2.0*1.534698e-18*slice_length/(gamma**3*beta**2))*(macroSize/bunch_length)/(bunch_radius**2)

slope_avg = 0.
slope2_avg = 0.
count = 0
for ip in range(b.getSize()):
	x = b.x(ip)
	y = b.y(ip)
	xp = b.xp(ip)	
	yp = b.yp(ip)	
	z = b.z(ip)
	dE = b.dE(ip)
	r = math.sqrt(x*x+y*y)
	p = math.sqrt(xp*xp+yp*yp)
	if( r > 0.5*bunch_radius and r < 0.9*bunch_radius):
		slope_avg += p/r
		slope2_avg += (p/r)**2
		count += 1 
slope_avg /= count
slope2_avg /= count
slope_err = math.sqrt((slope2_avg - slope_avg*slope_avg))
print "particles slope delta(p)/r            = %12.5g +- %12.5g "%(slope_avg,slope_err)
print "particles slope delta(p)/r from theory= %12.5g"%slope_theory 

nStep = 300
rho_arr = []
for ix in range(nStep+1):
	x = ix*bunch_radius/nStep
	y = 0.
	r = math.sqrt(x*x+y*y)
	rho = rhoGrid.getValue(x,y)
	rho_arr.append((r,rho))
	
phi_arr = []
phi_00 = phiGrid.getValue(0.,0.)
for ix in range(2*nStep+1):
	x = (ix-nStep)*bunch_radius/nStep
	y = 0.
	r = math.sqrt(x*x+y*y)
	phi = phiGrid.getValue(x,y) - phi_00
	phi_arr.append((x,phi))
	
long_arr = []
for ix in range(nStep+1):
	z = ix*bunch_length/nStep - bunch_length/2.0
	long_rho = longGrid.getValue(z)/longGrid.getStepZ()
	long_arr.append((z,long_rho))
	
long_grad_arr = []
for ix in range(nStep+1):
	z = ix*bunch_length/nStep - bunch_length/2.0
	long_grad_rho = longDerivGrid.getValue(z)/longGrid.getStepZ()
	long_grad_arr.append((z,long_grad_rho))
		

#-------------------------------------------------	
#this is the example of using the Gnuplot package
#-------------------------------------------------
import Gnuplot
gRho = Gnuplot.Gnuplot()
gRho.title('Rho SC')
gRho('set data style line')
gRho.plot(rho_arr)

gPhi = Gnuplot.Gnuplot()
gPhi.title('Phi SC')
gPhi('set data style line')
gPhi.plot(phi_arr)

gLong = Gnuplot.Gnuplot()
gLong.title('Long SC')
gLong('set data style line')
gLong.plot(long_arr)

gLongGrad = Gnuplot.Gnuplot()
gLongGrad.title('Long Gradient SC')
gLongGrad('set data style line')
gLongGrad.plot(long_grad_arr)

raw_input('Please press return to stop:\n')
#-------------------------------------------------	

sys.exit(1)

print "Stop."

