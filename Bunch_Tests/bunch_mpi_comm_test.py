import sys

from bunch import Bunch

import orbit_mpi
from orbit_mpi import mpi_comm

#-----------------------------------------------------
#get and set MPI Comm communicators for bunch
#-----------------------------------------------------

print "Start."

b = Bunch()

comm_local = b.getMPIComm()
print "comm name=",orbit_mpi.MPI_Comm_get_name(comm_local)

#create group and comm from wider comm and sub-group
group = orbit_mpi.MPI_Comm_group(comm_local)  
comm_new = orbit_mpi.MPI_Comm_create(mpi_comm.MPI_COMM_WORLD,group)
#set name to the comm
orbit_mpi.MPI_Comm_set_name(comm_new,"new comm")

b.setMPIComm(comm_new)
print "new comm name=",orbit_mpi.MPI_Comm_get_name(b.getMPIComm())

print "Stop."

#orbit_mpi.finalize("The test is done!")

# --- Memory leak test ---
b = Bunch()
count = 0
while(True):
	count = count + 1
	comm_local = b.getMPIComm()
	group = orbit_mpi.MPI_Comm_group(comm_local) 
	comm_new = orbit_mpi.MPI_Comm_create(mpi_comm.MPI_COMM_WORLD,group)
	orbit_mpi.MPI_Comm_set_name(comm_new,"new comm")
	b.setMPIComm(comm_new)
	s = orbit_mpi.MPI_Comm_get_name(b.getMPIComm())
	if(count % 10000 == 0):
		print "count=",count
