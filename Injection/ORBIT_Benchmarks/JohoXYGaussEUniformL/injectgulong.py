import math
import sys
from bunch import Bunch
from foil import Foil
from collimator import Collimator
from injection import InjectParts
from injection import JohoTransverse, GULongDist
from orbit.utils.orbit_mpi_utils import bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit
print "Start."

xmin = -0.050
xmax = 0.050
ymin = -0.050
ymax = 0.050

foilparams = (xmin, xmax, ymin, ymax)

#------------------------------
#Bunch init
#------------------------------
b = Bunch()
runName = "Test_Injection"

b.mass(0.93827231)
b.macroSize(1.0e+1)
energy = 1.0 #Gev
b.getSyncParticle().kinEnergy(energy)

lostfoilbunch = Bunch()
lostfoilbunch.addPartAttr("LostParticleAttributes") 

#------------------------------
#Initial Distribution Functions
#------------------------------
sp = b.getSyncParticle()

order = 3.
alphax = 0.063
betax = 10.209
alphay = 0.063
betay = 10.776
emitlim = 0.00152 * 2*(order + 1) * 1e-6
xcenterpos = 0.0468
xcentermom = 0.001
ycenterpos = 0.0492
ycentermom = -0.00006
tailfrac = 0
tailfac = 1
zlim = 120. * 248./360.
zmin = -zlim
zmax = zlim
emean = 1.0
esigma = 0.001
etrunc = 1
emin = 0.997
emax = 1.003

sp = b.getSyncParticle()
xFunc = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom, tailfrac, tailfac)
yFunc = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom, tailfrac, tailfac)
lFunc = GULongDist(zmin, zmax, sp, emean, esigma, etrunc, emin, emax)

#------------------------------
# Inject some particles
#------------------------------

nparts = 10000
inject = InjectParts(nparts, b, lostfoilbunch, foilparams, xFunc, yFunc, lFunc)
inject.addParticles()

bunch_pyorbit_to_orbit(248.0, b, "pybunch.dat")
lostfoilbunch.dumpBunch("pylostbunch.dat")
print "Stop."
quit()


