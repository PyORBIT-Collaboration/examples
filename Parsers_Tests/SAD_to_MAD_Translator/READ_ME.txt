There are python scripts and mad8 executable for preparation of the MAD
output files.
The rcs_sad_to_mad_translator.py will translate the lattice structure
from a SAD input file to the MAD format lattice file.
The rcs_.601.mad has the command to run mad and to prepare
TWISS and LATTICE. These files will be used by ORBIT and by the python
script to prepare TEAPOT-like imput for ORBIT.

The generation rcs_2.601.lat (MAD Lattice)
>./rcs_sad_to_mad_translator.py rcs_2.601.sad rcs_2.601.lat

Run mad8 to get all mad output files and pictures:
>./mad8.linux.bin < rcs_2.601_noda.mad > out
or
>./mad8.linux.bin < rcs_2.601.mad > out

The rcs_2.601_noda.mad does not include a tune matching,
and therefore it does not have the same optics as for rcs_2.601.mad.

The MR080221i.sad file is the Main Ring SAD file for JAERY. 


====== MULT elements in SAD =============
If L > 0 in SAD it will be a drift in MAD
