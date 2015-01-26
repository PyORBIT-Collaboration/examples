##############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice and add injection 
##############################################################
import sys
import pickle
import math
import numpy as np

from bunch import Bunch


# lattice, teapot class
from orbit.teapot import teapot
from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import DriftTEAPOT
from orbit.teapot import TEAPOT_MATRIX_Lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer
from bunch import Bunch, BunchTwissAnalysis

# add injection node, distribution generators 
from orbit.injection import TeapotInjectionNode
from orbit.injection import addTeapotInjectionNode

#from injection import InjectParts
from orbit.injection import UniformLongDist
from orbit.bunch_generators import TwissContainer, TwissAnalysis
from orbit.bunch_generators import KVDist1D, KVDist2D


# add apertures
from orbit.aperture import addTeapotApertureNode
from orbit.aperture import TeapotApertureNode, CircleApertureNode, EllipseApertureNode, RectangleApertureNode
from orbit.aperture import addCircleApertureSet, addEllipseApertureSet, addRectangleApertureSet
from aperture import Aperture



# trans space charge
from orbit.space_charge.sc2p5d import scAccNodes, scLatticeModifications
from orbit.space_charge.sc3d import scAccNodes, scLatticeModifications
from orbit.space_charge.directforce2p5d import directforceAccNodes, directforceLatticeModifications
from spacecharge import SpaceChargeForceCalc2p5D


from orbit_utils import BunchExtremaCalculator


import scipy as sp
from scipy import constants


    
#==================redif. of physic cont. to PATRIC convention========
clight = constants.c
eps0 = constants.value("electric constant")
qe = constants.e
mp = constants.m_p
me = constants.m_e
re = constants.value("classical electron radius")
mp_Mev = constants.value("proton mass energy equivalent in MeV")
rp = 1.53469e-18
#==================redif. of physic cont. to PATRIC convention========
    
    
#=====set up bunch stuff============
NPIC = 5000
bunch = Bunch()
A = 238
Z = 28
energy = 11.4e-3*A # Gev
intensity = 1.4e10
        
emittance_x = 500e-6
emittance_y = 500e-6

bunch.mass(0.93827231*A)
bunch.charge(Z)
bunch.macroSize(intensity/NPIC/A)
bunch.getSyncParticle().kinEnergy(energy)


sp = bunch.getSyncParticle()
gamma0 = sp.gamma()
beta0 = sp.beta()
E0 = A*mp*gamma0*pow(beta0*clight,2)/qe*1e-9


paramsDict = {}
lostbunch = Bunch()
paramsDict["lostbunch"]=lostbunch
paramsDict["bunch"]= bunch
lostbunch.addPartAttr("LostParticleAttributes") 

filename = "fodo_out.lat"
#filename = "fodo_simply_out.lat"
teapot_latt = teapot.TEAPOT_Lattice()
teapot_latt.readMADX(filename,"cella")

circum = teapot_latt.getLength()

for node in teapot_latt.getNodes():
    print "node=", node.getName()," s start,stop = %4.3f %4.3f "%teapot_latt.getNodePositionsDict()[node]
    if filename == "fodo_out.lat":
        if node.getName() == "qd":
            print "poles", node.getParam("poles"), "skews", node.getParam("skews"), "kls", node.getParam("kls"), "kq", node.getParam("kq")
        if node.getName() == "Aperture":
            print "aperture limitation", node.getParam("aperture")
        if node.getName() == "qf":
            print "poles", node.getParam("poles"), "skews", node.getParam("skews"), "kls", node.getParam("kls"), "kq", node.getParam("kq")


            
