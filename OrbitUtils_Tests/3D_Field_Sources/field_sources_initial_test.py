#------------------------------------------------------
# This is an example of Field Sources of electric and magnetic
# fields for 3D tracking using RK4_Tracker 
#-------------------------------------------------------

import sys
import math
import time

from orbit_utils import field_sources
from field_sources import QuadFieldSource
from field_sources import DipoleFieldSource
from field_sources import MagnetFieldSourceGrid3D

from orbit_utils import Matrix, PhaseVector
from spacecharge import Grid3D

transfCoordsMatrix = Matrix(4,4)
transfCoordsMatrix.unit()

#----------------------------------------------
#   Quad set get parameters test
#----------------------------------------------

quad_field_source = QuadFieldSource()
quad_field_source.length(1.5)
print "quad length =",quad_field_source.length()
quad_field_source.gradient(3.0)
print "quad gradient =",quad_field_source.gradient()

quad_field_source.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = quad_field_source.transormfMatrix()
transfCoordsMatrix.unit()

#----------------------------------------------
#   Dipole set get parameters test
#----------------------------------------------

dipole_field_source = DipoleFieldSource()
dipole_field_source.fieldsXYZ(0.1,0.2,0.3)
print "dipole fieldsXYZ = ",dipole_field_source.fieldsXYZ()
dipole_field_source.sizesXYZ(0.5,0.6,0.7)
print "dipole sizesXYZ = ",dipole_field_source.sizesXYZ()

dipole_field_source.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = dipole_field_source.transormfMatrix()
transfCoordsMatrix.unit()

#----------------------------------------------
#   Grid3D field source set get parameters test
#----------------------------------------------

bx_grid3d = Grid3D(10,10,10)
by_grid3d = Grid3D(10,10,10)
bz_grid3d = Grid3D(10,10,10)

magnet_field3D = MagnetFieldSourceGrid3D(bx_grid3d,by_grid3d,bz_grid3d)

magnet_field3D.transormfMatrix(transfCoordsMatrix)
transfCoordsMatrix = magnet_field3D.transormfMatrix()
transfCoordsMatrix.unit()

sys.exit(0)

#---------------------------------------------------
# Memory leak and performance test
#---------------------------------------------------

print "----- start loop --------"


start_time = time.clock()

count = 1
while( 1 < 2):

	quad_field_source = QuadFieldSource()
	dipole_field_source = DipoleFieldSource()
		
	transfCoordsMatrix = Matrix(4,4)
	transfCoordsMatrix.unit()
	dipole_field_source.transormfMatrix(transfCoordsMatrix)
	transfCoordsMatrix = dipole_field_source.transormfMatrix()
	transfCoordsMatrix.unit()	
	
	quad_field_source.transormfMatrix(transfCoordsMatrix)
	transfCoordsMatrix = quad_field_source.transormfMatrix()
	
	bx_grid3d = Grid3D(10,10,10)
	by_grid3d = Grid3D(10,10,10)
	bz_grid3d = Grid3D(10,10,10)
	magnet_field3D_1 = MagnetFieldSourceGrid3D(bx_grid3d,by_grid3d,bz_grid3d)
	magnet_field3D_2 = MagnetFieldSourceGrid3D(bx_grid3d,by_grid3d,bz_grid3d)
	magnet_field3D_3 = MagnetFieldSourceGrid3D(bx_grid3d,by_grid3d,bz_grid3d)
	magnet_field3D_1.transormfMatrix(transfCoordsMatrix)
	transfCoordsMatrix = magnet_field3D_1.transormfMatrix()
	
	(Ex,Ey,Ez,Bx,By,Bz) = magnet_field3D_1.getFields(0.,0.,0.)
	(Ex,Ey,Ez,Bx,By,Bz) = quad_field_source.getFields(0.,0.,0.)
	(Ex,Ey,Ez,Bx,By,Bz) = dipole_field_source.getFields(0.,0.,0.)
	
	if(count % 10000 == 0):
		stop_time = time.clock() - start_time
		performance = stop_time/count
		print "count =",count,"   performance loop/sec= %10.3g "%performance
		
	count += 1


print "Stop"