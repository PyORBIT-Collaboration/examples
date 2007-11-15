import orbit_mpi
from orbit_mpi import mpi_comm

#------------------------------------------------------
# This example includes folowing MPI functions:
# orbit_mpi.MPI_Comm_create(comm_old,group)
# orbit_mpi.MPI_Comm_group(comm)
# orbit_mpi.MPI_Comm_dup(comm)
# orbit_mpi.MPI_Comm_set_name(comm,name)
# orbit_mpi.MPI_Comm_get_name((comm)
# orbit_mpi.MPI_Comm_compare(comm1,comm2)
#
# From the mpi_comm sub-package: 
# mpi_comm.MPI_Comm()
# 
# comm = mpi_comm.MPI_COMM_WORLD
# comm = mpi_comm.MPI_COMM_SELF
# comm = mpi_comm.MPI_COMM_NULL
# 
#------------------------------------------------------

rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)


comm_world = mpi_comm.MPI_COMM_WORLD
comm_self  = mpi_comm.MPI_COMM_SELF
comm_null  = mpi_comm.MPI_COMM_NULL

comm = mpi_comm.MPI_Comm()

print "MPI_COMM_WORLD=",comm_world
print "MPI_COMM_SELF=",comm_self
print "MPI_COMM_NULL=",comm_null
print "comm=",comm

if(orbit_mpi.MPI_Comm_compare(comm,comm_world) == orbit_mpi.MPI_IDENT):
	if(rank == 0): print "new commumocator is IDENT to MPI_COMM_WORLD"

#copy operation
comm_copy = orbit_mpi.MPI_Comm_dup(comm)

#create group
group = orbit_mpi.MPI_Comm_group(comm)

#create comm from wider comm and sub-group  
comm_new = orbit_mpi.MPI_Comm_create(mpi_comm.MPI_COMM_WORLD,group)

#set name to the comm
orbit_mpi.MPI_Comm_set_name(comm_new,"new comm")
print "name=",orbit_mpi.MPI_Comm_get_name(comm_new)

if(rank == 0): print "=========================================="

orbit_mpi.finalize("The test is done!")

#--------------------------
# CHECK memory leak
#--------------------------
count = 0
while(1 < 2):
	count += 1
	comm_copy = orbit_mpi.MPI_Comm_dup(comm)
	group = orbit_mpi.MPI_Comm_group(comm_copy)
	comm_new = orbit_mpi.MPI_Comm_create(mpi_comm.MPI_COMM_WORLD,group)
	orbit_mpi.MPI_Comm_set_name(comm_new,"new comm")
	name = orbit_mpi.MPI_Comm_get_name(comm_new)
	if(count % 10000 == 0 and rank == 0):
		print "i=",count
