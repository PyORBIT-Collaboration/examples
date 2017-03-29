This directory includes two methods for calculating Transport Matrices
between two bunches. 

1. The first method is presented in the file 
transport_matr_gen_test_1.py

The method is based on the correspondence between two bunches according to the 
Particles Id Attribute. In the beginning we assign Id numbers to each particle 
in the bunch, create copy of this bunch,and track the particles in the bunch 
through the known transport matrix. Then we compare particles in the copy of 
the original bunch to the particles in the bunch after the tracking. The result is
a transport matrix. The particles could be lost during the tracking.

2. The example of the second method is in 
transport_matr_gen_test_2.py

Here we save the initial coordinates in the Init Coordinates Attributes of 
the bunch, and after tracking compare the new coordinates with the initial
ones.