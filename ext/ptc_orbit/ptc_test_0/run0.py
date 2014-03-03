import sys
import os
import math

from orbit.teapot_base import TPB

from orbit.utils import orbitFinalize

from orbit.lattice import AccLattice, AccNode,\
     AccActionsContainer, AccNodeBunchTracker

from bunch import Bunch
from orbit.utils.orbit_mpi_utils import\
     bunch_orbit_to_pyorbit, bunch_pyorbit_to_orbit

from orbit.teapot import TEAPOT_Lattice
from orbit.teapot import BaseTEAPOT

from libptc_orbit import *

from ext.ptc_orbit import PTC_Lattice
from ext.ptc_orbit import PTC_Node
from ext.ptc_orbit.ptc_orbit import setBunchParamsPTC, readAccelTablePTC,\
     readScriptPTC, updateParamsPTC, synchronousSetPTC, synchronousAfterPTC,\
     trackBunchThroughLatticePTC, trackBunchInRangePTC

PTC_File = "ptc_data_test_0.txt"

length_of_name = len(PTC_File)
ptc_init_(PTC_File, length_of_name - 1)

Lattice = PTC_Lattice("TestCase")

Lattice.readPTC(PTC_File)

print Lattice.getLength(), Lattice.betax0, Lattice.betay0, Lattice.alphax0, Lattice.alphay0, Lattice.etax0, Lattice.etapx0
PhaseLength = Lattice.getLength()
print Lattice.getLength(), PhaseLength

"""
for node in Lattice.getNodes():
    print node.getType(), node.getLength(), node.getParam("node_index"), node.getParam("betax"), node.getParam("betay"), node.getParam("alphax"), node.getParam("alphay"), node.getParam("etax"), node.getParam("etapx")
"""

b = Bunch()
print "Read Bunch."
runName = "PTC Test"

setBunchParamsPTC(b)
kin_Energy = b.getSyncParticle().kinEnergy()
print kin_Energy, b.charge(), b.mass()

total_macroSize=1.0e+10

bunch_orbit_to_pyorbit(Lattice.getLength(), kin_Energy, "bunch_ini.dat", b)

nParticlesGlobal = b.getSizeGlobal()
b.macroSize(total_macroSize/nParticlesGlobal)
print nParticlesGlobal, b.macroSize()

updateParamsPTC(Lattice, b)

synchronousSetPTC(-1)

"""
synchronousAfterPTC(-1)
"""

Turns = 1
for i in range(Turns):
    print i
    trackBunchThroughLatticePTC(Lattice, b, PhaseLength)

"""
trackBunchInRangePTC(Lattice, b, PhaseLength, 0, 500)

trackBunchInRangePTC(Lattice, b, PhaseLength, 501, 932)

trackBunchInRangePTC(Lattice, b, PhaseLength, 0, 10)

trackBunchInRangePTC(Lattice, b, PhaseLength, 11, 20)

trackBunchInRangePTC(Lattice, b, PhaseLength, 21, 932)
"""

b.dumpBunch("bunch_temp.dat")
bunch_pyorbit_to_orbit(Lattice.getLength(), b, "bunch_final.dat")

readScriptPTC("ptc_data_test_0.txt")
