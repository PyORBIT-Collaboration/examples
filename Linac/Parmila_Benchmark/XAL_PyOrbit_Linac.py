#   author:     T. Gorlov 
#   date:       08/24/2012
#   
#   
#   Run the script with jython:
#
#>  jython XAL_PyOrbit_Linac.py
#
#   First, it will run the jython part of the script to generate HEBT_lattice.txt data file and bunch0.txt
#   After that it will authomatically run the PyOrbit part of the script for optics calculation 
#
#   The emittance parameters of input bunch bunch0.txt have XAL definition of emittance 
#   Pay attention that lattice.readMAD("input.mad",'HEBT1_HEBT2') gives incorrect result  
#   because magnetic field parameters have sign uncertainty at this time
#
#   output files:
#   PyORBIT_rms.txt xal_rms.txt     are output files with the following structure:
#   z sigma_x sigma_y sigma_z Tkin


import sys

if (len(sys.argv) == 1):
    
    from gov.sns.xal.smf.data import *
    from gov.sns.ca import *
    from gov.sns.xal.smf import *
    from java.io import *
    from gov.sns.xal.smf import *
    from gov.sns.xal.smf import *
    from java.util import *
    from gov.sns.tools.plot import *
    from gov.sns.xal.smf.impl import *
    from gov.sns.ca import *

    from gov.sns.xal.model.probe import *
    from gov.sns.xal.model.scenario import *
    from gov.sns.xal.model.probe import *
    from gov.sns.xal.model.probe.traj import *
    from gov.sns.xal.model.alg import *
    from gov.sns.xal.model import *
    from gov.sns.xal.model.pvlogger import *
    from gov.sns.tools.beam import Twiss
    import sys, os, math, random
    

    acc  = XMLDataManager.loadDefaultAccelerator()

    lst = ArrayList()
    sclm = acc.getSequence("SCLMed")
    sclh = acc.getSequence("SCLHigh")
    hebt1 = acc.getSequence("HEBT1")
    hebt2 = acc.getSequence("HEBT2")

    lst.add(sclm)
    lst.add(sclh)
    lst.add(hebt1)
    lst.add(hebt2)

    seq = AcceleratorSeqCombo("LINAC", lst)

    model = Scenario.newScenarioFor(seq)
    model.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)


    probe = EnvelopeProbe()
    T = 0.1864805061222586
    #T = 0.925                      #(GeV)
    m = 0.938256 + 2*0.000511                     #(GeV)

    E = m + T
    pz = math.sqrt(E*E - m*m)
    q = -1
    beta = pz/E
    g = E/m



    Tx = Twiss(-1.42, 9.60, 0.298e-6)
    Ty = Twiss(-0.265, 6.23, 0.189e-6)
    Tz = Twiss(0.565, 5.33, 0.928e-6)
    
    probe.initFromTwiss([Tx,Ty,Tz])
    probe.setBeamCurrent(0)
    probe.setSpeciesCharge(q)
    probe.setSpeciesRestEnergy(m*1.0e9)
    probe.setKineticEnergy(T*1.0e9)
    probe.setAlgorithm(EnvelopeTracker())

    model.setProbe(probe)
    model.resync()
    model.run()

    trajectory = model.getTrajectory()

    f1 = open("xal_rms.txt","w")  
    for i in range(trajectory.numStates()):
        state = trajectory.stateWithIndex(i)       
        z = state.getPosition()
        s = state.phaseCorrelation()
        T = state.getKineticEnergy()

        sx = s.getSigmaX()
        sy = s.getSigmaY()
        sz = s.getSigmaZ()


        f1.write(str(z)+"\t"+str(sx)+"\t"+str(sy)+"\t"+str(sz)+"\t"+str(T)+"\n")

    f1.close()



    #-----------------------------------bunch generate----------------------------------------

    b = open("bunch_0.txt",'w')    
    b.write(("\
# PARTICLE_ATTRIBUTES_CONTROLLERS_NAMES     \n\
# BUNCH_ATTRIBUTE_DOUBLE charge   %i    \n\
# BUNCH_ATTRIBUTE_DOUBLE classical_radius   1.5347e-18  \n\
# BUNCH_ATTRIBUTE_DOUBLE macro_size   0     \n\
# BUNCH_ATTRIBUTE_DOUBLE mass   %f    \n\
#  SYNC_PART_COORDS 0 0 0  x, y, z positions in [m] \n\
#  SYNC_PART_MOMENTUM 0 0 %f  px, py, pz momentum component in GeV/cn   \n\
#  SYNC_PART_X_AXIS 1 0 0  nxx, nxy, pxz - x-axis ort coordinatesn  \n\
#  info only: energy of the synchronous particle [GeV] = %f  \n\
#  info only: momentum of the synchronous particle [GeV/c] = %f     \n\
#  info only: beta=v/c of the synchronous particle = %f   \n\
#  info only: gamma=1/sqrt(1-(v/c)**2) of the synchronous particle = %f    \n\
#  SYNC_PART_TIME 0  time in [sec]  \n\
#  x[m] px[rad] y[m] py[rad] z[m]  (pz or dE [GeV])\n"%(q, m, pz, T, pz, beta, g)).replace('#', '%'))

    def getCoords(alpha, beta, emit_rms, u_cutoff):
        """ Returns u and u' [m] and [rad] """

        while 1:
            u = random.gauss(0.,math.sqrt(beta*emit_rms))
            up = random.gauss(0.,math.sqrt(emit_rms/beta))
            if(u*u+beta*beta*up*up < u_cutoff*u_cutoff):
                up = up - alpha*u/beta
                return u,up

    for i in range(100000):
    
        x,xp = getCoords(Tx.getAlpha(), Tx.getBeta(), Tx.getEmittance(), 0.05)
        y,yp = getCoords(Ty.getAlpha(), Ty.getBeta(), Ty.getEmittance(), 0.05)
        z,zp = getCoords(Tz.getAlpha(), Tz.getBeta(), Tz.getEmittance(), 0.05) 

        b.write(str(x)+" "+str(xp)+" "+str(y)+" "+str(yp)+" "+str(z)+" "+str(g*g*beta*pz*zp)+"\n")


    b.close()
    
    os.system('${ORBIT_ROOT}/bin/pyORBIT XAL_PyOrbit_Linac.py calc')
    sys.exit()    
      
    
if (sys.argv[1] == "calc"):
    print "calculate"

    import sys
    import math

    from orbit.sns_linac import SimplifiedLinacParser
    from orbit.sns_linac import LinacLatticeFactory, LinacAccLattice

    from bunch import Bunch

    from orbit.lattice import AccLattice, AccNode, AccActionsContainer

    parser = SimplifiedLinacParser("../SNS_Linac_XML/sns_linac_hebt.xml")
    linacTree = parser.getLinacStructureTree()
    print "======================================="
    print "Total length=",linacTree.getLength()
    print "======================================="
    sequences = linacTree.getSeqs()
    totalLength = 0.
    for seq in sequences:
            totalLength +=  seq.getLength()	
            print "seq=",seq.getName()," L=",seq.getLength(),"  total length=",totalLength    
    
            
    lattFactory = LinacLatticeFactory(linacTree)
    accLattice = lattFactory.getLinacAccLattice(["SCLMed","SCLHigh","HEBT1", "HEBT2"])
    #accLattice = lattFactory.getLinacAccLattice(["UserSequence",])
    #accLattice = lattFactory.getLinacAccLattice(["SCLMed"])
    
    b = Bunch()
    b.readBunch("bunch_0.txt")

    paramsDict = {"test_pos":0.,"count":0}
    actionContainer = AccActionsContainer("Test Design Bunch Tracking")
    
    def action_entrancee(paramsDict):
            node = paramsDict["node"]
            length = node.getLength()
            pos = paramsDict["test_pos"] + length
            paramsDict["test_pos"] = pos	
            bunch = paramsDict["bunch"]

            eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3	
            if(node.getName().find(":Rg") >= 0):
                    paramsDict["count"]	+= 1
                    s = " %5d     %25s     %4.5f     %5.3f  "%(paramsDict["count"],node.getName(),(pos - length/2),eKin)
                    #outF.write(s+"\n")
                    print s
    
    
    print "Start."

    file_out = open("pyORBIT_rms.txt","w")

    def action_entrance(paramsDict):
            if(isinstance(paramsDict["parentNode"],AccLattice)):

                node = paramsDict["node"]
                pos = paramsDict["test_pos"]
                bunch = paramsDict["bunch"]

                eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3



                def get_rms():
                    x0 = 0
                    y0 = 0
                    z0 = 0

                    for i in range(bunch.getSize()):
                        x0 += bunch.x(i)
                        y0 += bunch.y(i)
                        z0 += bunch.y(i)

                    x0 = x0/bunch.getSize()
                    y0 = y0/bunch.getSize()
                    z0 = z0/bunch.getSize()

                    sx = 0
                    sy = 0
                    sz = 0

                    for i in range(bunch.getSize()):
                        sx += (x0 - bunch.x(i))**2
                        sy += (y0 - bunch.y(i))**2
                        sz += (z0 - bunch.z(i))**2

                    sx = math.sqrt(sx/bunch.getSize())
                    sy = math.sqrt(sy/bunch.getSize())
                    sz = math.sqrt(sz/bunch.getSize())

                    return sx,sy,sz

                sx, sy, sz = get_rms()

                file_out.write(str(pos)+"\t"+str(sx)+"\t"+str(sy)+"\t"+str(sz)+"\t"+str(eKin) +"\n")

                length = node.getLength()
                pos = paramsDict["test_pos"] + length
                paramsDict["test_pos"] = pos


    
    actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)      
    
    accLattice.trackDesignBunch(b, paramsDict = paramsDict)

    accLattice.trackBunch(b, paramsDict = paramsDict, actionContainer = actionContainer)
    
    #accLattice.trackBunch(b, paramsDict = paramsDict)

    #b.dumpBunch("bunch_ls_out.txt")   


    file_out.close()

    sys.exit() 