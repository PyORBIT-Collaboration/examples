#!/bin/bash

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the python script>"
    exit $E_BADARGS
fi

pyORBIT $1 $2 $3 $4 $5 $6 $7 $8 $9 
