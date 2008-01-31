import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_group

#------------------------------------------------------
# This example includes folowing MPI functions:
# orbit_mpi.MPI_Comm_group(comm)
# orbit_mpi.MPI_Group_size(group)
# orbit_mpi.MPI_Group_rank(group)
# orbit_mpi.MPI_Group_union(group1,group2)
#
#
# There are more these functions!
#
# From the mpi_group sub-package: 
# mpi_group.MPI_Group()
# mpi_group.MPI_GROUP_NULL
# mpi_group.MPI_GROUP_EMPTY
# 
# There is NO MPI_GROUP_WORLD !!!
#------------------------------------------------------
group_null = mpi_group.MPI_GROUP_NULL
group_empty = mpi_group.MPI_GROUP_EMPTY

group_world = orbit_mpi.MPI_Comm_group(orbit_mpi.mpi_comm.MPI_COMM_WORLD)

size = orbit_mpi.MPI_Group_size(group_world)
rank = orbit_mpi.MPI_Group_rank(group_world)
print "rank=",rank," size=",size

group = orbit_mpi.MPI_Group_union(group_empty,group_world)

if(rank == 0): print "=========================================="

#orbit_mpi.finalize("The test is done!")

#--------------------------
# CHECK memory leak
#--------------------------
count = 0
while(1 < 2):
	count += 1
	group_empty = mpi_group.MPI_GROUP_EMPTY
	group = orbit_mpi.MPI_Group_union(group_empty,group_world)
	if(count % 10000 == 0 and rank == 0):
		print "i=",count
