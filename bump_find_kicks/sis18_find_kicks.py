##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice and add injection 
##############################################################
import sys
import math
import numpy as np
from pylab import *
from matplotlib.font_manager import FontProperties


from bunch import Bunch

# lattice, teapot class
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import TEAPOT_MATRIX_Lattice


from orbit.bumps import simpleBump, close_orbit_bumps
from orbit.orbit_correction import orbit, correction

from orbit.matrix_lattice import PlotBeamline


#=====set up bunch stuff============
bunch = Bunch()
print "Read Bunch."

A = 1
Z = 1
energy = 11.4e-3*A # Gev
intensity = 1.4e10
emittance_x = 50e-6
emittance_y = 50e-6
bunch.mass(0.93827231)
bunch.getSyncParticle().kinEnergy(energy)
#=====set up bunch stuff============

	
print "Generate Lattice."
lattice = TEAPOT_Lattice("no_sc_lattice")

lattice.readMAD("sis18_inj.lat","SIS18_MID") # lattice start at injection point plus half circumference 
	


# get the lattice function for the kick and bpms, jet only horizontal 
find = close_orbit_bumps()

xc0 = 70e-3
xcs0 = 7.0e-3 


# set variable, list [[node.name], [variable, start value], ...]
kicker1 = [["S11MB1"], ["kx",-0.00111]]
kicker2 = [["S12MB2"], ["kx",0.00]]
kicker3 = [["S01MB3"], ["kx",0.0333]]
kicker4 = [["S03MB4"], ["kx", -0.00111]]
variables = [kicker1, kicker2, kicker3, kicker4]

# set constraints, list [[node.name], [coordinate (x or xp), value], ...]
con1 = [["S12DX5H"], ["x", xc0]]
con2 = [["S11MB1"], ["x", 0.0]]
con3 = [["S03MB4"], ["x", 0.0]]
con4 = [["S12DX5H"], ["xp", xcs0]]
constraints = [con1, con2, con3, con4]


find.lattice_function(lattice, bunch, variables, constraints)

kick_leastsq = find.find_kicks(method="LeastSq")

kick = kick_leastsq
nodes = lattice.getNodes()
for node in nodes:
	if node.getName() == "S11MB1":
		node.setParam("kx",kick[0])
	if node.getName() == "S12MB2":
		node.setParam("kx",kick[1])
	if node.getName() == "S01MB3":
		node.setParam("kx",kick[2])
	if node.getName() == "S03MB4":
		node.setParam("kx",kick[3])


OrbitX, OrbitY = orbit(lattice,bunch).get_orbit()
x = []
y = []
s = []
for i in xrange(len(OrbitX)):
	s.append(OrbitX[i][0])
	x.append(OrbitX[i][1]*1e3)
	y.append(OrbitX[i][2]*1e3)
	
	
# matplotlib parameter
rc('text', usetex=False) 
rc('font', family='serif')
rc('font', size=18)
rc('axes', linewidth=2)
rc('lines', linewidth=3)
rc('figure', figsize=(8,6))



Pline = PlotBeamline(lattice)

min_s = 72
max_s = 165

figure(figsize=[10, 10])
ax0 = plt.subplot2grid((9, 1), (0, 0))
ax1 = plt.subplot2grid((9, 1), (1, 0), rowspan=4)
ax2 = plt.subplot2grid((9, 1), (5, 0), rowspan=4)

title="SEC11-SEC03"

Pline.plot_beamline(ax0,l1=min_s,l2=max_s, title=title, septum_name="S12DX5H")

ax1.plot(s,x,'r-', label="x")
ax2.plot(s,y,'b-', label="y")
#setp(ax1.get_xticklabels(), visible=False)

ax1.set_ylabel(r'$x$ [mm]')
ax2.set_ylabel(r'$x_p$ [mrad]')
ax2.set_xlabel('s [m]')


ax1.set_xlim(min_s,max_s)
ax2.set_xlim(min_s,max_s)
ax2.set_ylim(-20,20)

fontP = FontProperties()
fontP.set_size(35)

subplots_adjust(bottom=0.08,left=0.15,hspace=0.5,right=0.95, top=0.98)

#savefig("plot.pdf")

show()
exit()