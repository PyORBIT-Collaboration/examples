import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

#------------------------------------------------------
# This example includes folowing MPI functions:
# orbit_mpi.MPI_Barrier(comm)
# orbit_mpi.MPI_Send(data,data_type,dest_rank,tag,comm)
# orbit_mpi.MPI_Recv(data_type,source,tag,comm) => data_recv
# orbit_mpi.MPI_Bcast(data,data_type,main_rank,comm) => data_res
# orbit_mpi.MPI_Allreduce(data,data_type,op,comm) => data_res
# 
#------------------------------------------------------

comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD
rank = orbit_mpi.MPI_Comm_rank(comm)
size = orbit_mpi.MPI_Comm_size(comm)

if(size < 2):
	message = "Number of CPUs should be more than 2!"
	message = message + " nCPUs=" + str(size)
	orbit_mpi.finalize(message)

#data could be scalar values also
arr = (1,2,3)
arr = "TEST"
arr = (1.1,2.1,3.1)

data_type = mpi_datatype.MPI_INT
data_type = mpi_datatype.MPI_CHAR
data_type = mpi_datatype.MPI_DOUBLE

if(rank == 0):
	orbit_mpi.MPI_Send(arr,data_type,1,111,comm)

if(rank == 1):
	arr_rcv = orbit_mpi.MPI_Recv(data_type,0,111,comm)
	print "rank=",rank," RECV arr=",arr_rcv


orbit_mpi.MPI_Barrier(comm)
if(rank == 0): print "=========================================="

op = mpi_op.MPI_SUM
arr_new = orbit_mpi.MPI_Allreduce(arr,data_type,op,comm)
if(rank == 0): print "allreduce arr=",arr_new

main_rank = 0
arr_new = orbit_mpi.MPI_Bcast(arr,data_type,main_rank,comm)
if(rank == 1): print "bcast arr=",arr_new

#orbit_mpi.finalize("The test is done!")

#--------------------------
# CHECK memory leak
#--------------------------
count = 0
while(1 < 2):
	count += 1
	if(rank == 0):
		orbit_mpi.MPI_Send(arr,data_type,1,111,comm)
	if(rank == 1):
		arr_rcv = orbit_mpi.MPI_Recv(data_type,0,111,comm)
	arr_new = orbit_mpi.MPI_Allreduce(arr,data_type,op,comm)
	arr_new = orbit_mpi.MPI_Bcast(arr,data_type,main_rank,comm)
	if(count % 1000 == 0 and rank == 0):
		print "i=",count
