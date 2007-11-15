import orbit_mpi
import sys

cpu  = orbit_mpi.MPI_Get_processor_name()
rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
tm = orbit_mpi.MPI_Wtime()
tk = orbit_mpi.MPI_Wtick()
print "rank=",rank," size=",size, " name=", cpu," time=",tm," tick=",tk

comm_world = orbit_mpi.mpi_comm.MPI_COMM_WORLD
comm = orbit_mpi.mpi_comm.MPI_Comm()
orbit_mpi.MPI_Comm_dup(comm_world,comm)
orbit_mpi.MPI_Comm_set_name(comm,"my_communicator")
print "comm_world=",orbit_mpi.mpi_comm.MPI_COMM_WORLD
print "comm_self=",orbit_mpi.mpi_comm.MPI_COMM_SELF
print "comm=",comm
name = orbit_mpi.MPI_Comm_get_name(comm)
print "name comm=",name

group1 = orbit_mpi.mpi_group.MPI_Group()
group2 = orbit_mpi.mpi_group.MPI_Group()
orbit_mpi.MPI_Group_incl(group1,[-5.1,2,3],group2)

sys.exit(1)

count = 0
while(1 < 2):
	count += 1
	orbit_mpi.MPI_Comm_set_name(comm,"my_communicator"+str(count))
	name = orbit_mpi.MPI_Comm_get_name(comm)
	if(count % 10000 == 0):
		print "i=",count

#orbit_mpi.finalize("Error Message!")
orbit_mpi.finalize()
