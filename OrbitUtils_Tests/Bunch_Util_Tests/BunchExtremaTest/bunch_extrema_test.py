#-----------------------------------------------------
# Calculates the min and max coordinates for the bunch
#-----------------------------------------------------

import sys
import math

from orbit_utils import BunchExtremaCalculator

import orbit_mpi
from orbit_mpi import mpi_comm

from bunch import Bunch

bunch_extrema_cal = BunchExtremaCalculator()

rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(orbit_mpi.mpi_comm.MPI_COMM_WORLD)

b = Bunch()

nParts = 3
for i in xrange(nParts):
	b.addParticle(0.1+i+rank,0.2+i+rank,0.3+i+rank,0.4+i+rank,0.5+i+rank,0.6+i+rank)
	

(xMin,xMax,yMin,yMax,zMin,zMax) = bunch_extrema_cal.extremaXYZ(b)

if(rank == 0):
	print "xMin xMax = ",xMin," ",xMax
	print "yMin yMax = ",yMin," ",yMax
	print "zMin zMax = ",zMin," ",zMax

