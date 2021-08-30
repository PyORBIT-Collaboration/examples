#==================================
# This script is a test for the matrix determinant 
# calculations
#==================================


import sys
import time

from orbit_utils import Matrix

def printM(m):
	print "----matrix--- size=",m.size()
	for i in xrange(m.size()[0]):
		for j in xrange(m.size()[1]):
			print ("m(" + str(i) + "," + str(j)+")="+str(m.get(i,j)) + " "),
		print ""	

def setM(arr):
	n = len(arr)
	m = Matrix(n,n)
	for i in xrange(n):
		for j in xrange(n):
			m.set(i,j,arr[i][j])
	return m
		
print "Start."

"""

arr = [[],[]]
arr[0] = [ 1., 0.]
arr[1] = [ 0., 3.]

arr = [[],[],[]]
arr[0] = [ 1., 0., 0.]
arr[1] = [ 0., 2., 0.]
arr[2] = [ 0., 0., 3.]

arr = [[],[],[],[]]
arr[0] = [ 1., 0., 0., 0.]
arr[1] = [ 0., 2., 0., 0.]
arr[2] = [ 0., 0., 3., 0.]
arr[3] = [ 0., 0., 0., 4.]

"""

arr = [[],[],[],[],[]]
arr[0] = [ 1., 0., 0., 0., 0.]
arr[1] = [ 0., 2., 0., 0., 0.]
arr[2] = [ 0., 0., 3., 0., 0.]
arr[3] = [ 0., 0., 0., 4., 0.]
arr[4] = [ 0., 0., 0., 0., 5.]         
             
mtrx = setM(arr)
print "==================================="
printM(mtrx)
print "==================================="
print "Det(mtrx)=",mtrx.det()

sys.exit(0)

count = 1
start_time = time.time()
while(1<2):
	det = mtrx.det()
	if(count % 1000000 == 0):
		print "count =",count," speed [calc/sec] = ",count/(time.time()-start_time)
	count += 1
	
