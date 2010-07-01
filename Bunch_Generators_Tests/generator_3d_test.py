#!/usr/bin/env python

#--------------------------------------------------------
# The script will test the gerenrators from bunch_generators
#--------------------------------------------------------

import math
import random
import sys

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

n = 10000
#---------------------------------------------
# KV 3D  
#---------------------------------------------
twissX = TwissContainer(alpha = 1., beta = 2., emittance = 3.)
twissY = TwissContainer(alpha = 2., beta = 3., emittance = 4.)
twissZ = TwissContainer(alpha = 3., beta = 4., emittance = 5.)
dist = KVDist3D(twissX,twissY,twissZ)
twiss_analysis = TwissAnalysis(3)
for i in range(n):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp,z,zp))
	
print "================================================="
print "KV 3D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "......................................................................"
print "Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g "%twissZ.getAlphaBetaGammaEmitt()
print "Generated Z  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(2)
print "================================================="

#---------------------------------------------
# Water Bag 3D 
#---------------------------------------------
dist = WaterBagDist3D(twissX,twissY,twissZ)
twiss_analysis.init()
for i in range(n):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp,z,zp))
	
print "================================================="
print "Water Bag 3D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "......................................................................"
print "Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g "%twissZ.getAlphaBetaGammaEmitt()
print "Generated Z  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(2)
print "================================================="

#---------------------------------------------
# Gauss 3D 
#---------------------------------------------
dist = GaussDist3D(twissX,twissY,twissZ)
twiss_analysis.init()
for i in range(n):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp,z,zp))
	
print "================================================="
print "Gauss 2D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "......................................................................"
print "Twiss     Z  %12.5g  %12.5g   %12.5g    %12.5g "%twissZ.getAlphaBetaGammaEmitt()
print "Generated Z  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(2)
print "================================================="


