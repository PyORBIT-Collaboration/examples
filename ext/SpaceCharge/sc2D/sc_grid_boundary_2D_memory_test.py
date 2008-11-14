#-----------------------------------------------------
#Creates Grid2D and Boundary2D to check the memory leak
#-----------------------------------------------------
import sys

from spacecharge import Grid2D
from spacecharge import Boundary2D

print "Start."

count = 0
while( 1 < 2):
	grid = Grid2D(100,200)
	boundary = Boundary2D(100,200,10.,20.,64,"Ellipse",20)
	grid.setBoundary2D(boundary)
	grid = Grid2D(boundary)
	grid.setBoundary2D(boundary)
	grid = Grid2D(-5.0,5.0,100,-10.,10.,200)
	grid.setBoundary2D(boundary)
	count = count + 1
	if(count % 100 == 0): print "count=",count

print "Stop."

sys.exit(1)
