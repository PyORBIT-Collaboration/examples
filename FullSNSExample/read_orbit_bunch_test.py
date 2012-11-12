##############################################################
# This script reads the ORBIT particles file and distribute it 
# among CPUs 
##############################################################

import time
import math

from bunch import Bunch
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

import orbit_mpi
from orbit_mpi import mpi_op
from orbit_mpi import mpi_datatype

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)	

b = Bunch()

#
ORBIT_file_name = "Bm_Parts_1_0"
ORBIT_NEW_file_name = "new_orbit_file.dat"
pyORBIT_DUMP_file_name = "pyorbit_bunch_dump_test.dat"

kineticEnergy = 1.0
ringLength = 248.0

time_start = time.clock()

bunch_orbit_to_pyorbit(ringLength, kineticEnergy, ORBIT_file_name,b)
print "rank=",rank," n_parts=",b.getSize()
bunch_pyorbit_to_orbit(ringLength, b, ORBIT_NEW_file_name)
time_exec = time.clock() - time_start
if(rank == 0): print "time[sec]=",time_exec

b.dumpBunch(pyORBIT_DUMP_file_name)

x_sum = 0.
xp_sum = 0.
y_sum = 0.
yp_sum = 0.
phi_sum = 0.
dE_sum = 0.	

for i in range(b.getSize()):
	x = b.x(i)
	xp = b.xp(i)
	y = b.y(i)
	yp = b.yp(i)
	z = b.z(i)
	dE = b.dE(i)	
	x_sum += x
	xp_sum += xp
	y_sum += y
	yp_sum += yp
	phi_sum += z
	dE_sum += dE	
		
var_arr = (x_sum,xp_sum,y_sum,yp_sum,phi_sum,dE_sum)
var_arr = orbit_mpi.MPI_Allreduce(var_arr,mpi_datatype.MPI_DOUBLE,mpi_op.MPI_SUM,comm)
(x_sum,xp_sum,y_sum,yp_sum,phi_sum,dE_sum) = var_arr
if(rank == 0):
	print "================ parallel ============="
	print "x_sum   =",x_sum
	print "xp_sum  =",xp_sum
	print "y_sum   =",y_sum
	print "yp_sum  =",yp_sum
	print "phi_sum =",phi_sum
	print "dE_sum  =",dE_sum
		
if(rank == 0):
	x_sum = 0.
	xp_sum = 0.
	y_sum = 0.
	yp_sum = 0.
	phi_sum = 0.
	dE_sum = 0.	
	file_in = open(ORBIT_file_name,"r")
	ln = file_in.readline().strip()
	
	while ln:
		res_arr = ln.strip().split()
		x  = float(res_arr[0])/1000.
		xp = float(res_arr[1])/1000.
		y  = float(res_arr[2])/1000.
		yp = float(res_arr[3])/1000.
		z  = float(res_arr[4])*ringLength/(2*math.pi)
		dE = float(res_arr[5])
		x_sum += x
		xp_sum += xp
		y_sum += y
		yp_sum += yp
		phi_sum += z
		dE_sum += dE		
		ln = file_in.readline().strip()
		
	print "===========single CPU ================"
	print "x_sum   =",x_sum
	print "xp_sum  =",xp_sum
	print "y_sum   =",y_sum
	print "yp_sum  =",yp_sum
	print "phi_sum =",phi_sum
	print "dE_sum  =",dE_sum
	file_in.close()

#----------------------------------------------------------
#  read new pyORBIT file and analyze
#----------------------------------------------------------
b = Bunch()
b.readBunch(pyORBIT_DUMP_file_name)

x_sum = 0.
xp_sum = 0.
y_sum = 0.
yp_sum = 0.
phi_sum = 0.
dE_sum = 0.	

for i in range(b.getSize()):
	x = b.x(i)
	xp = b.xp(i)
	y = b.y(i)
	yp = b.yp(i)
	z = b.z(i)
	dE = b.dE(i)	
	x_sum += x
	xp_sum += xp
	y_sum += y
	yp_sum += yp
	phi_sum += z
	dE_sum += dE	
		
var_arr = (x_sum,xp_sum,y_sum,yp_sum,phi_sum,dE_sum)
var_arr = orbit_mpi.MPI_Allreduce(var_arr,mpi_datatype.MPI_DOUBLE,mpi_op.MPI_SUM,comm)
(x_sum,xp_sum,y_sum,yp_sum,phi_sum,dE_sum) = var_arr
if(rank == 0):
	print "================ read new pyORBIT file and analyze ============="
	print "x_sum   =",x_sum
	print "xp_sum  =",xp_sum
	print "y_sum   =",y_sum
	print "yp_sum  =",yp_sum
	print "phi_sum =",phi_sum
	print "dE_sum  =",dE_sum

