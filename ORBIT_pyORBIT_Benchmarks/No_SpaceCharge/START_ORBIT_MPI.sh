#!/bin/bash

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N-CPUs>"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N CPUs>"
    exit $E_BADARGS
fi

mpirun -np $2 $ORBIT_PATH/$CPU/ORBIT_MPI $1
