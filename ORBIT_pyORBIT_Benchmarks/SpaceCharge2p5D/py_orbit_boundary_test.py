
from spacecharge import Boundary2D

# boundary 
nBoundaryPoints = 34
N_FreeSpaceModes = 10
boundary_radius = 0.11
boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Circle",2.0*boundary_radius)

print "Circle is done."


nBoundaryPoints = 34
N_FreeSpaceModes = 10
x_dim = 2.
y_dim = 2.

boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Ellipse",x_dim,y_dim)

print "Ellipse is done."

nBoundaryPoints = 1200
N_FreeSpaceModes = 256
x_dim = 2.
y_dim = 2.

boundary = Boundary2D(nBoundaryPoints,N_FreeSpaceModes,"Rectangle",x_dim,y_dim)

print "Rectangle is done."