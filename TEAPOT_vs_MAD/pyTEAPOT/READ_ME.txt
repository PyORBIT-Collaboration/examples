Then scripts in this directory:

mad_matrix_reading.py - reads the MAD output file and print the one
                        turn matrix. Users can print all elements'
                        matrices one by one with the names of elements.

teapot_stand_alone_node_test.py - the test for generating the transport matrices
                        for stand alone TEAPOT elements. Users can check that the
                        resulting matrices are what they should be.

teapot_matrix_calculation.py - different approaches for transport matrix 
                        calculations. 
                        1. An action in a action container and a track method of 
                           TEAPOT nodes
                        2. A transport matrix generation for the whole TEAPOT lattice
                        3. A trackBunch method for each TEAPOT node
												

matrix_lattice_test.py - MATRIX_Lattice class is presented. The example for for 
                        calculation of ring parameters and twiss functions for 
                        the ring and accelerator lines are presented 
  

matrix_lattice_test.out - the results of running of the matrix_lattice_test.py script to compare.
