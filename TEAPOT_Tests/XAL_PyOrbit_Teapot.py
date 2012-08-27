#   author:     T. Gorlov 
#   date:       08/24/2012
#   
#   
#   Run the script with jython:
#
#>  jython teapot_HEBT.py
#
#   First, it will run the jython part of the script to generate HEBT_lattice.txt data file and bunch0.txt
#   After that it will authomatically run the PyOrbit part of the script for optics calculation 
#
#   The emittance parameters of input bunch bunch0.txt have XAL definition of emittance 
#   Pay attention that lattice.readMAD("input.mad",'HEBT1_HEBT2') gives incorrect result  
#   because magnetic field parameters have sign uncertainty at this time
#
#   output files:
#   PyOrbit.txt xal_rms.txt     are output files with the following structure:
#   z sigma_x sigma_y sigma_z Disp_x Disp_dx


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

    initial_data = PVLoggerDataSource(15990732)

    lst = ArrayList()
    hbt1 = acc.getSequence("HEBT1")
    hbt2 = acc.getSequence("HEBT2")

    lst.add(hbt1)
    lst.add(hbt2)

    seq = AcceleratorSeqCombo("HEBT", lst)

    model = Scenario.newScenarioFor(seq)
    model.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)

    model_pv = initial_data.setModelSource(seq, model)

    probe = EnvelopeProbe()
    
    T = 0.925                       #(GeV)
    m = 0.938256 + 0.000511         #(GeV)

    E = m + T
    pz = math.sqrt(E*E - m*m)
    q = -1
    beta = pz/E
    g = E/m

    Tx = Twiss(2.533267, 7.55, 731.6004e-9)
    Ty = Twiss(-3.673985, 21.16748, 522.4011e-9)
    Tz = Twiss(-115.54787, 7098.4685, 2.557899e-9)
    
    probe.initFromTwiss([Tx,Ty,Tz])
    probe.setBeamCurrent(0)
    probe.setSpeciesCharge(q)
    probe.setSpeciesRestEnergy(m*1.0e9)
    probe.setKineticEnergy(T*1.0e9)
    probe.setAlgorithm(EnvelopeTracker())

    model_pv.setProbe(probe)
    model_pv.resync()   
    model_pv.run()

    trajectory = model_pv.getTrajectory()

    f1 = open("xal_rms","w")  
    for i in range(trajectory.numStates()):
        state = trajectory.stateWithIndex(i)       
        z = state.getPosition()
        s = state.phaseCorrelation()

        sx = s.getSigmaX()
        sy = s.getSigmaY()
        sz = s.getSigmaZ()
        Dx = (s.getElem(0, 5) - s.getElem(0, 6)*s.getElem(5, 6))/(s.getElem(5, 5)-s.getElem(5, 6)*s.getElem(5, 6))/(g*g);
        Dpx = (s.getElem(1, 5) - s.getElem(1, 6)*s.getElem(5, 6))/(s.getElem(5, 5)-s.getElem(5, 6)*s.getElem(5, 6))/(g*g);

        f1.write(str(z)+"\t"+str(sx)+"\t"+str(sy)+"\t"+str(sz)+"\t"+str(Dx)+"\t"+str(Dpx)+"\n")

    f1.close()



    
    
    lattice = model_pv.getLattice()


    f2 = open("HEBT_lattice.txt",'w',0)
    for i in range(lattice.getChild(0).getChildCount()):
        node = lattice.getChild(0).getChild(i)

        f2.write(node.getType() + "\t" + node.getId() + "\t")

        if(node.getType() == "IdealMagQuad"):
            f2.write(str(node.getLength()) + "\t" + str(node.getMagField()))
        if(node.getType() == "IdealMagWedgeDipole2"):
            sb = node.getChild(1)
            f2.write(str(node.getLength()) + "\t" +str(-sb.getFieldIndex()*sb.compDesignCurvature()*sb.compDesignCurvature()) +"\t" + str(node.getDesignBendingAngle()) + "\t" + str(node.getEntrPoleAngle()) + "\t" + str(node.getExitPoleAngle()))    
        if(node.getType() == "IdealDrift"):
            f2.write(str(node.getLength()))
        if(node.getType() == "IdealMagSteeringDipole"):
            f2.write(str(node.getMagField()))



        f2.write("\n")

    f2.close()

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
    
    os.system('${ORBIT_ROOT}/bin/pyORBIT XAL_PyOrbit_Teapot.py calc')
    sys.exit()    
      
    
if (sys.argv[1] == "calc"):
    print "calculate"

    from orbit.teapot import teapot
    from orbit.lattice import AccLattice, AccActionsContainer
    from bunch import Bunch
    import math



    lattice = teapot.TEAPOT_Lattice()


    b = Bunch()
    b.readBunch("bunch_0.txt")
    pz = b.getSyncParticle().momentum()
    q = b.charge()
    m = b.mass()
    E = math.sqrt(pz*pz + m*m)
    k = E/(pz*pz)


    #lattice.readMAD("input.mad",'HEBT1_HEBT2')

    f = open("HEBT_lattice.txt",'r')
    while 1:
        arr = f.readline().split()
        if (len(arr) == 0): break
        elem = None

        if (arr[0] == "IdealMagSteeringDipole"):
            elem = teapot.KickTEAPOT(arr[1])   
            if("DCH" in arr[1]):
                elem.addParam("kx",0)
            if("DCV" in arr[1]):
                elem.addParam("ky",0)

        if (arr[0] == "IdealDrift"):
            elem = teapot.DriftTEAPOT(arr[1])
            elem.setLength(float(arr[2]))

        if (arr[0] == "IdealMagQuad"):
            elem = teapot.QuadTEAPOT(arr[1])
            elem.setLength(float(arr[2]))        
            elem.addParam("kq",q*0.299792458*float(arr[3])/pz)

        if (arr[0] == "IdealMagWedgeDipole2"):
            elem = teapot.BendTEAPOT(arr[1]) 
            elem.setLength(float(arr[2]))
            elem.addParam("kls",[-q*float(arr[3])*elem.getLength()])     
            elem.addParam("theta",float(arr[4]))
            elem.addParam("ea1",float(arr[5]))
            elem.addParam("ea2",float(arr[6]))        
            elem.addParam("poles",[1])
            elem.addParam("skews",[0])

        if (arr[0] == "Marker"):
            elem = teapot.NodeTEAPOT(arr[1])      



        lattice.addNode(elem)

    lattice.initialize()


    print "Start."

    file_out = open("pyORBIT_rms.txt","w")


    def action_entrance(paramsDict):
            if(isinstance(paramsDict["parentNode"],AccLattice)):

                node = paramsDict["node"]
                pos = paramsDict["test_pos"]
                bunch = paramsDict["bunch"]


                if (node.getName() == "LS_IP_2"):
                    bunch.dumpBunch("bunch_ls.txt")


                def get_rms():
                    x0 = 0
                    y0 = 0
                    z0 = 0
                    Ep0 = 0
                    xp0 = 0


                    for i in range(bunch.getSize()):
                        x0 += bunch.x(i)
                        y0 += bunch.y(i)
                        z0 += bunch.y(i)
                        Ep0 += bunch.pz(i)
                        xp0 += bunch.px(i)


                    x0 = x0/bunch.getSize()
                    y0 = y0/bunch.getSize()
                    z0 = z0/bunch.getSize()
                    Ep0 = Ep0/bunch.getSize()
                    xp0 = xp0/bunch.getSize()

                    sx = 0
                    sy = 0
                    sz = 0
                    Ep2 = 0
                    xEp = 0
                    xpEp = 0


                    for i in range(bunch.getSize()):
                        sx += (x0 - bunch.x(i))**2
                        sy += (y0 - bunch.y(i))**2
                        sz += (z0 - bunch.z(i))**2
                        Ep2 += (Ep0 - bunch.pz(i))**2
                        xEp += (x0 - bunch.x(i))*(Ep0 - bunch.pz(i))
                        xpEp += (xp0 - bunch.px(i))*(Ep0 - bunch.pz(i))


                    sx = math.sqrt(sx/bunch.getSize())
                    sy = math.sqrt(sy/bunch.getSize())
                    sz = math.sqrt(sz/bunch.getSize())
                    Ep2 = Ep2/bunch.getSize()
                    xEp = xEp/bunch.getSize()
                    xpEp = xpEp/bunch.getSize()


                    sD = xEp/(Ep2*k)
                    sDp = xpEp/(Ep2*k)

                    return sx,sy,sz,sD,sDp

                sx, sy, sz, sD, sDp  = get_rms()



                file_out.write(str(pos)+"\t"+str(sx)+"\t"+str(sy)+"\t"+str(sz)+"\t"+str(sD)+"\t"+str(sDp) +"\n")

                #file_out.write(str(pos)+"\t"+str((bunch.x(1) - bunch.x(0))/0.000001)+"\n")
                length = node.getLength()
                pos = paramsDict["test_pos"] + length
                paramsDict["test_pos"] = pos




    actionContainer = AccActionsContainer("Test Design Bunch Tracking")                
    actionContainer.addAction(action_entrance, AccActionsContainer.ENTRANCE)	
    paramsDict = {"test_pos":0.,"count":0}
    lattice.trackBunch(b, paramsDict = paramsDict, actionContainer = actionContainer)

    file_out.close()

    print "=========================================="

    print "lattice length=",lattice.getLength()
    print "beta=",b.getSyncParticle().beta()
    print "TEAPOT time[sec]=",b.getSyncParticle().time()
    print "SIMPLE time[sec]=",lattice.getLength()/(b.getSyncParticle().beta()*2.99792458e+8)
    print "Stop."

    sys.exit()
    
