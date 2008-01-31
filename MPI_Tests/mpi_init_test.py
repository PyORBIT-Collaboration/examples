import orbit_mpi

#------------------------------------------------------
# This example includes folowing MPI functions
# MPI_Initialized()
# MPI_Get_processor_name()
# MPI_Comm_rank( MPI_Comm )
# MPI_Comm_size( MPI_Comm )
# MPI_Wtime()
# MPI_Wtick()
# finalize([message])
#------------------------------------------------------
# and MPI constants:
# MPI_IDENT
# MPI_CONGRUENT
# MPI_SIMILAR
# MPI_UNEQUAL
# MPI_UNDEFINED
# MPI_UNDEFINED_RANK
# MPI_SUCCESS
# MPI_ANY_SOURCE
# MPI_ANY_TAG
#------------------------------------------------------

mpi_init = orbit_mpi.MPI_Initialized()

cpu  = orbit_mpi.MPI_Get_processor_name()
rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
size = orbit_mpi.MPI_Comm_size(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
t = orbit_mpi.MPI_Wtime()
tick = orbit_mpi.MPI_Wtick()

#---------MPI Constant -----------------------
if(rank == 0):
	print "init=",mpi_init," rank=",rank," size=",size, " name=", cpu," time=",t," tick=",tick
	print "MPI_IDENT=",orbit_mpi.MPI_IDENT
	print "MPI_CONGRUENT=",orbit_mpi.MPI_CONGRUENT
	print "MPI_SIMILAR=",orbit_mpi.MPI_SIMILAR
	print "MPI_UNEQUAL=",orbit_mpi.MPI_UNEQUAL
	print "MPI_UNDEFINED=",orbit_mpi.MPI_UNDEFINED
	print "MPI_UNDEFINED_RANK=",orbit_mpi.MPI_UNDEFINED_RANK
	print "MPI_SUCCESS=",orbit_mpi.MPI_SUCCESS
	print "MPI_ANY_SOURCE=",orbit_mpi.MPI_ANY_SOURCE
	print "MPI_ANY_TAG=",orbit_mpi.MPI_ANY_TAG
else:
	print "init=",mpi_init," rank=",rank," size=",size, " name=", cpu

#orbit_mpi.finalize("Test Error Message!")
#orbit_mpi.finalize()

#--------------------------
# CHECK memory leak
#--------------------------
count = 0
while(1 < 2):
	count += 1
	mpi_init = orbit_mpi.MPI_Initialized()
	cpu  = orbit_mpi.MPI_Get_processor_name()
	rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
	size = orbit_mpi.MPI_Comm_size(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
	if(count % 10000 == 0 and rank == 0):
		print "i=",count
