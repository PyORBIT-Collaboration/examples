This is directory with benchmark between ORBIT_MPI and pyORBIT for tracking 
macro-particles through a simple FODO lattice.

All scripts here are NON-PARALLEL! Do not use them with more than 1 CPU!

Fist run pyORBIT script:
>./START.sh py_orbit_no_sc.py 1

It will create the orbit_mpi_bunch_input.dat file (which is input for ORBIT_MPI script ORBIT_MPI_NO_SC.sc) 
and orbit_mpi_bunch_from_pyorbut_output.dat file which is results from pyORBIT tracking translated to
ORBIT_MPI bunch format.

Then run ORBIT_MPI script ORBIT_MPI_NO_SC.sc:
>./START_ORBIT_MPI.sh ORBIT_MPI_NO_SC.sc 1

The result will be a file orbit_mpi_bunch_output.dat which should be compared with orbit_mpi_bunch_from_pyorbut_output.dat.

Today results are in files:
orbit_mpi_bunch_from_pyorbut_output_arch.dat
orbit_mpi_bunch_output_arh.dat
orbit_mpi_bunch_input_arch.dat

You can compare them with yours.

