#!/usr/bin/env python

#--------------------------------------------------------
# The script will test the gerenrators from bunch_generators
#--------------------------------------------------------

import math
import random
import sys

from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist2D, GaussDist2D, KVDist2D


n = 10000
#---------------------------------------------
# KV 2D  
#---------------------------------------------
twissX = TwissContainer(alpha = 1., beta = 2., emittance = 3.)
twissY = TwissContainer(alpha = 2., beta = 3., emittance = 4.)
dist = KVDist2D(twissX,twissY)
twiss_analysis = TwissAnalysis(2)
for i in range(n):
	(x,xp,y,yp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp))
	
print "================================================="
print "KV 2D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "================================================="


#---------------------------------------------
# Water Bag 2D  
#---------------------------------------------
dist = WaterBagDist2D(twissX,twissY)
twiss_analysis.init()
for i in range(n):
	(x,xp,y,yp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp))
	
print "================================================="
print "Water Bag 2D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "================================================="


#---------------------------------------------
# Gauss 2D  
#---------------------------------------------
dist = WaterBagDist2D(twissX,twissY)
twiss_analysis.init()
for i in range(n):
	(x,xp,y,yp) = dist.getCoordinates()
	twiss_analysis.account((x,xp,y,yp))
	
print "================================================="
print "Gauss 2D - done!"
print "                  alpha       beta [cm/rad]    gamma     emitt[cm*rad] "
print "Twiss     X  %12.5g  %12.5g   %12.5g    %12.5g "%twissX.getAlphaBetaGammaEmitt()
print "Generated X  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(0)
print "......................................................................"
print "Twiss     Y  %12.5g  %12.5g   %12.5g    %12.5g "%twissY.getAlphaBetaGammaEmitt()
print "Generated Y  %12.5g  %12.5g   %12.5g    %12.5g "%twiss_analysis.getTwiss(1)
print "================================================="


