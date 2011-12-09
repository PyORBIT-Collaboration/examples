#!/usr/bin/env python

#--------------------------------------------------------
# The script will test the gerenrators from bunch_generators
#--------------------------------------------------------

import math
import random
import sys

from orbit.bunch_generators import TwissContainer, KVDist1D, WaterBagDist1D, GaussDist1D
from orbit.bunch_generators import TwissAnalysis

# number of particles generated
n = 10000

#---------------------------------------------
# KV 1D  
#---------------------------------------------
twiss = TwissContainer(1.,2.,3.)
dist = KVDist1D(twiss)
twiss_analysis = TwissAnalysis(1)
for i in range(n):
	(x,y) = dist.getCoordinates()
	twiss_analysis.account((x,y))

print "================================================="
print "KV 1D - done!"
print "                  alpha       beta [m/rad]    gamma     emitt[m*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "================================================="

#---------------------------------------------
# Waterbag 1D  
#---------------------------------------------
dist = WaterBagDist1D(twiss)
twiss_analysis.init()
for i in range(n):
	(x,y) = dist.getCoordinates()
	twiss_analysis.account((x,y))

print "================================================="
print "Water bag 1D - done!"
print "                  alpha       beta [m/rad]    gamma     emitt[m*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "================================================="

#---------------------------------------------
# Gauss 1D  
#---------------------------------------------
dist = GaussDist1D(twiss, cut_off = 5.0)
twiss_analysis.init()
for i in range(n):
	(x,y) = dist.getCoordinates()
	twiss_analysis.account((x,y))

print "================================================="
print "Gauss 1D - done!"
print "                  alpha       beta [m/rad]    gamma     emitt[m*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "================================================="



