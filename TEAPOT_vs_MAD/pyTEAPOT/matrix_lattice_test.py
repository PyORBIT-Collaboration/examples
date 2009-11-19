##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and then it builds the MATRIX lattice
# with the first order transport matrices. It can calculate the ring
# parameters such as:
# tunes X,Y for the ring
# momentum compaction, gamma transition, slip factor
# twiss parameters along the accelerator line and for the ring
# chromaticities X,Y
# 
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import BaseTEAPOT
from orbit.teapot import RingRFTEAPOT
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.matrix_lattice import BaseMATRIX
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit_utils import Matrix
from bunch import Bunch

#---PRINT Function for Matrix
def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("(%1d,%1d)=% 6.5e "%(i,j,m.get(i,j))),
		print ""	

print "Start."

teapot_latt = teapot.TEAPOT_Lattice()
print "Read MAD."
teapot_latt.readMAD("../MAD/LATTICE","RING")
ring_length = teapot_latt.getLength()
print "Lattice=",teapot_latt.getName()," ring length [m] =",ring_length

Tkin = 1.0
bunch = Bunch()
bunch.getSyncParticle().kinEnergy(Tkin)

matrix_lattice = TEAPOT_MATRIX_Lattice(teapot_latt,bunch)
print "Lattice=",matrix_lattice.getName()," matrix lattice length [m] =",matrix_lattice.getLength()

one_turn_matrix = matrix_lattice.getOneTurnMatrix()
print "=============one turn pyORBIT matrix for MATRIX lattice===="
printM(one_turn_matrix)

print "=== ring parameters ==="
ring_par_dict = matrix_lattice.getRingParametersDict()
Tkin = ring_par_dict["Ekin [GeV]"]
T = ring_par_dict["period [sec]"]
freq = ring_par_dict["frequency [Hz]"]
eta = ring_par_dict["eta"]
gamma_trans = ring_par_dict["gamma transition"]
Tkin_trans = ring_par_dict["transition energy [GeV]"]
alpha_p = ring_par_dict["momentum compaction"] 
f_tune_x = ring_par_dict["fractional tune x"]
f_tune_y = ring_par_dict["fractional tune x"]

beta_x = ring_par_dict["beta x [m]"]
beta_y = ring_par_dict["beta y [m]"]
alpha_x = ring_par_dict["alpha x"]
alpha_y = ring_par_dict["alpha y"]
disp_x = ring_par_dict["dispersion x [m]"]
disp_y = ring_par_dict["dispersion y [m]"]
disp_prime_x = ring_par_dict["dispersion prime x"]
disp_prime_y = ring_par_dict["dispersion prime y"]

print " Ekin [GeV]= %5.3f   Period T [sec]= %12.5g   Freq. [Hz]= %12.5g "%(Tkin,T,freq)
print " Slip Factor eta= %12.5g  Gamma Trans.= %12.5g   Tkin Trans.= %12.5g "%(eta,gamma_trans,Tkin_trans)
print " Momentum compuction alpha_p= %12.5g "%(alpha_p)
print " Fractional tune x and y = %12.5g  %12.5g  "%(f_tune_x,f_tune_y)
print " twiss beta [m] x and y = %12.5g  %12.5g  "%(beta_x,beta_y)
print " twiss alpha x and y = %12.5g  %12.5g  "%(alpha_x,alpha_y)
print " dispersion [m] x and y = %12.5g  %12.5g  "%(disp_x,disp_y)
print " dispersion prime x and y = %12.5g  %12.5g  "%(disp_prime_x,disp_prime_y)

(tuneX, arrPosAlphaX, arrPosBetaX) = matrix_lattice.getRingTwissDataX()
(tuneY, arrPosAlphaY, arrPosBetaY) = matrix_lattice.getRingTwissDataY()
print "Tune x =",tuneX," graph points=",len(arrPosAlphaX)
print "Tune y =",tuneY," graph points=",len(arrPosAlphaY)

(arrDispX,arrDispPrimeX) = matrix_lattice.getRingDispersionDataX()
(arrDispY,arrDispPrimeY) = matrix_lattice.getRingDispersionDataY()

max_disp_x = 0.
max_disp_y = 0.
for (pos,disp) in arrDispX:
	if(abs(disp) > max_disp_x): max_disp_x = abs(disp)
for (pos,disp) in arrDispY:
	if(abs(disp) > max_disp_y): max_disp_y = abs(disp)	
print "max Disp. X [m] =",max_disp_x
print "max Disp. Y [m] =",max_disp_y

(chromX,chromY) = matrix_lattice.getChromaticitiesXY()
print "chromaticity X=",chromX
print "chromaticity Y=",chromY
"""
#-------------------------------------------------	
#this is the example of using the Gnuplot package
import Gnuplot
gBX = Gnuplot.Gnuplot()
gBX.title('beta x [m]')
gBX('set data style line')
gBX.plot(arrPosBetaX)

gBY = Gnuplot.Gnuplot()
gBY.title('beta y [m]')
gBY('set data style line')
gBY.plot(arrPosBetaY)

gDX = Gnuplot.Gnuplot()
gDX.title('dispersion x [m]')
gDX('set data style line')
gDX.plot(arrDispX)

#gDY = Gnuplot.Gnuplot()
#gDY.title('dispersion y [m]')
#gDY('set data style line')
#gDY.plot(arrDispY)

raw_input('Please press return to stop:\n')
#-------------------------------------------------	
"""
print "Stop."
sys.exit(1)


