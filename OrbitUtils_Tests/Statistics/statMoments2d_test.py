import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_utils import StatMoments2D
#----------------------------------------------------------------
# This test uses StatMoments2D to calculate statistical moments 
# and the Twiss parameters assuming we have 2D phase space
#----------------------------------------------------------------


rank = orbit_mpi.MPI_Comm_rank(orbit_mpi.mpi_comm.MPI_COMM_WORLD)
mpi_size = orbit_mpi.MPI_Comm_size(orbit_mpi.mpi_comm.MPI_COMM_WORLD)

if(rank == 1): print "MPI size=",mpi_size

statMoment2D = StatMoments2D()

statMoment2D.account(1.,2.)
statMoment2D.account(2.,3.)
statMoment2D.account(3.,1.)

#----if it is parallel the synchronization will sum all CPUs in communicator
#    by default communicator is comm_world
#statMoment2D.synchronizeMPI()
statMoment2D.synchronizeMPI(mpi_comm.MPI_COMM_WORLD)

if(rank == 1):
	
	max_order = statMoment2D.getMaxOrder()
	print "maxOrder = ",max_order
	for ui in range(max_order + 1):
		for upi in range(max_order + 1):
			print "moment( %2d , %2d )="%(ui,upi),statMoment2D.getStatMoment(ui,upi)
	
	print "count =",statMoment2D.getCount()
	
	print "emittance =",statMoment2D.getEmittance()
	print "Twiss alpha =",statMoment2D.getAlpha()
	print "Twiss beta  =",statMoment2D.getBeta()
	print "Twiss gamma =",statMoment2D.getGamma()
	print "Twiss rms(u)  =",statMoment2D.getRmsU()	
	print "Twiss rms(up) =",statMoment2D.getRmsUP()
