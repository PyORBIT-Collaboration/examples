#!/usr/bin/env python

#--------------------------------------------------------
# The script will test the gerenrators from bunch_generators
#--------------------------------------------------------

import math
import random
import sys

from bunch import Bunch
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import WaterBagDist3D, GaussDist3D, KVDist3D

n = 300000
#---------------------------------------------
# Water Bag 3D
#---------------------------------------------
twissX = TwissContainer(alpha = 0.045, beta = 12.355, emittance = 30e-6)
twissY = TwissContainer(alpha = 0.046, beta = 12.044, emittance = 30e-6)
twissZ = TwissContainer(alpha = 0, beta = 100000., emittance = .01)

#dist = WaterBagDist3D(twissX,twissY,twissZ)
dist = KVDist3D(twissX,twissY,twissZ)
#dist = GaussDist3D(twissX,twissY,twissZ)

twiss_analysis = TwissAnalysis(3)
twiss_analysis.init()
file_out = open("KV.dat","w")
for i in range(n):
	(x,xp,y,yp,z,zp) = dist.getCoordinates()
	file_out.write(str(x) + " " + str(xp) + " " + str(y) + " " + str(yp) + " "+ str(z) + " " + str(zp) + "\n")
	twiss_analysis.account((x,xp,y,yp,z,zp))	
file_out.close()
