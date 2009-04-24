##############################################################
# This script will read the matrices from the MAD output file.
# The one turn matrix will be printed at the end.
# Of course, you have to run MAD first.
##############################################################

import math
import sys

from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from orbit_utils import Matrix

#---PRINT Function for Matrix
def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("(%1d,%1d)=% 6.5e "%(i,j,m.get(i,j))),
		print ""	
		
print "Start."

inF = open("../MAD/mad_output")

s = inF.readline()

s_arr = []
stack_size = 8
elem_count = 0

lattice_arr = []

while ( s != ""):
	if(len(s_arr) < stack_size):
		s_arr.append(s)
	else:
		#======= we have the full stack of lines ===========
 		if(s_arr[1].find("Element transfer matrix:") >= 0):
			res_arr = s_arr[0].split()
			name = res_arr[len(res_arr) - 2][1:]
			elem_count = elem_count + 1
			#print "i=",elem_count," name=",name
			matrix = Matrix(6,6)
			for i in range(6):
				res_arr = s_arr[i+2].split()
				for j in range(6):
					matrix.set(i,j, float(res_arr[j]))
			#printM(matrix)
			lattice_arr.append((name,matrix))
		for it in range(stack_size-1):
			s_arr[it] = s_arr[it+1]
		s_arr[stack_size-1] = s	
	#print "debug s=",s
	s = inF.readline()

inF.close()

mad_matrix = Matrix(6,6)
mad_matrix.unit()
n_max = len(lattice_arr)
#n_max = 74
for i in range(n_max):
	(name,matrix) = lattice_arr[i]
	#print "i=",i," name =",name	
	#printM(mad_matrix)
	mad_matrix = matrix.mult(mad_matrix)	
	#printM(matrix)
	#printM(mad_matrix)


print "=========One turn MAD matrix======== MAD n_elements=",n_max 
printM(mad_matrix)

print "Stop."


