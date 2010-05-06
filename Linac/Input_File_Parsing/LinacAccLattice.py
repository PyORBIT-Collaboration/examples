"""
The general Linac Accelerator lattice. It tracks the bunch through the linac 
accelerator nodes. We cannot use here TEAPOT nodes from the TEAPOT package 
directly because we use the field as a parameter for quads and dipole correctors
instead of k1 = 1/(B*rho)*(dB/dr). The RF Cavities are different from the ring RF.
"""
import os
import math

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize, NamedObject, ParamsDictObject

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import the linac structure tree with all sequences and nodes, but without drifts
from linac_parser import LinacStructTree

class LinacAccLattice(AccLattice):
	"""
	The subclass of the AccLattice class. In the beginning the lattcie is empty.
	"""
	def __init__(self, name = None):
		AccLattice.__init__(self,name)
		self.__rfCavities = []
		self.__sequences = []
		
	def initialize(self):
		"""
		Method. Initializes the linac lattice, child node structures, and calculates 
		the one turn matrix.
		"""
		AccLattice.initialize(self)	
		
	
	def getSubLattice(self, index_start = -1, index_stop = -1,):
		"""
		It returns the new LinacAccLattice with children with indexes 
		between index_start and index_stop inclusive. 
		What about seqences and RF cavities ?????
		"""
		return self._getSubLattice(LinacAccLattice(),index_start,index_stop)

	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the lattice.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		paramsDict["bunch"] = bunch
		bunch.getSyncParticle().time(0.)
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)

	def trackDesignBunch(self, bunch_in):
		"""
		This will track the design bunch through the linac and set up RF Cavities times of
		arrivals.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Bunch Tracking")
		bunch = Bunch()
		bunch_in.copyEmptyBunchTo(bunch)
		bunch.getSyncParticle().time(0.)	
		paramsDict["bunch"] = bunch
		
		def trackDesignBunch(paramsDict):
			node = paramsDict["node"]
			node.trackDesignBunch(paramsDict)
			
		actionContainer.addAction(trackDesignBunch, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(trackDesignBunch, AccActionsContainer.BODY)

	def addRF_Cavity(self,cav):
		if(isinstance(cav, RF_Cavity) == True):
			self.__rfCavities.append(cav)
		else:
			msg = "The LinacAccLattice, method addRF_Cavity(cav)!"
			msg = msg + os.linesep
			msg = msg + "cav is not a subclass of RF_Cavity."
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)			

	def getRF_Cavity(self,name):
		""" Returns the cavity instance according to the name """
		for cav in self.__rfCavities:
			if(name == cav.getName()):
				return cav
		return None
		
	def getRF_Cavities(self):
		""" Returns the array with RF cavities. """
		return self.__rfCavities

	def addSequence(self,seq):
		if(isinstance(cav, Sequence) == True):
			self.__sequences.append(seq)
		else:
			msg = "The LinacAccLattice, method addSequence(seq)!"
			msg = msg + os.linesep
			msg = msg + "seq is not a subclass of Sequence."
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)			

	def getSequence(self,name):
		""" Returns the sequence instance according to the name """
		for seq in self.__sequences:
			if(name == seq.getName()):
				return seq
		return None
		
	def getSequences(self):
		""" Returns the array with sequences. """
		return self.__sequences

#--------------------------------------------------------------------------
#     Linac Lattice Factory
#--------------------------------------------------------------------------
class LinacLatticeFactory():
	""" . """
	def __init__(self, ltree):
		if(isinstance(ltree, LinacStructTree) != True):
			msg = "The LinacLatticeFactory constructor: you have to specify the LinacStructTree instance as input!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)	
		self.ltree = ltree
		#We need to compare positions, lengths etc. This is our delta
		self.zeroDistance = 0.00005
		
	def getLinacAccLattice(self,names):
		"""
		Returns the linac accelerator lattice for specified sequence names.
		"""
		if(len(names) < 1):
			msg = "The LinacLatticeFactory method getLinacAccLattice(names): you have to specify the names array!"
			msg = msg + os.linesep
			msg = msg + "Stop."
			msg = msg + os.linesep
			orbitFinalize(msg)				
		seqencesLocal = self.ltree.getSeqs()
		seqencesLocalNames = []
		for seq in seqencesLocal:
			seqencesLocalNames.append(seq.getName())
		ind_start = seqencesLocalNames.index(names[0])
		seqNames = seqencesLocalNames[ind_start:ind_start+len(names)]
		for name in seqNames:
			if(name not in names):
				msg = "The LinacLatticeFactory method getLinacAccLattice(names): sequence names array is wrong!"
				msg = msg + os.linesep
				msg = msg + "existing names=" + str(seqencesLocalNames)
				msg = msg + os.linesep
				msg = msg + "sequence names="+str(names)
				orbitFinalize(msg)
		sequences = self.ltree.getSeqs()[ind_start:ind_start+len(names)]
		#----make linac lattice
		linacAccLattice = LinacAccLattice(self.ltree.getName())
		#There are the folowing possible types of elements in the linac tree:
		#QUAD - quadrupole
		#RFGAP - RF Gap
		#DCH - horizontal dipole corrector
		#DCV - vertical dipole corrector
		#----------------------------------------------------------------------
		# The DRIFTS will be generated additionally and put into right places
		#----------------------------------------------------------------------
		accNodes = []
		accSeqs = []
		accRF_Cavs = [] 
		seqPosition = 0.
		for seq in sequences:
			accSeq = Sequence(seq.getName())
			accSeq.setLength(float(seq.getLength()))
			accSeq.setPosition(seqPosition)
			seqPosition = seqPosition + accSeq.getLength()
			accSeqs.append(accSeq)
			nodes = seq.getNodes()
			rf_cav_names = []
			thinNodes = []
			for node in nodes:
				node.setParam("pos",float(node.getParam("pos")))
				#------------QUAD-----------------
				if(node.getType() == "QUAD"):
					accNode = Quad(node.getName())
					accNode.setLength(float(node.getLength()))
					accNode.setParamsDict(node.getParamsDict())
					accSeq.addNode(accNode)
				#------------RF_Gap-----------------	
				elif(node.getType() == "RFGAP"):
					accNode = BaseRF_Gap(node.getName())
					accNode.setParamsDict(node.getParamsDict())
					accNode.setLength(float(node.getParam("gapLength")))
					accNode.setParam("length",accNode.getLength())
					accSeq.addNode(accNode)
					accNode.setParam("amp",float(node.getParam("amp")))
					accNode.setParam("E0TL",float(node.getParam("E0TL")))
					rf_cav_name = node.getParam("parent")
					if(rf_cav_name not in rf_cav_names):
						rf_cav_names.append(rf_cav_name)
						accRF_Cav = RF_Cavity(rf_cav_name)
						accRF_Cavs.append(accRF_Cav)
						accRF_Cav.setDesignPhase(accNode.getParam("firstPhase"))
						accRF_Cav.setPhase(accNode.getParam("firstPhase"))
						accRF_Cav.setDesignAmp(1.)
						accRF_Cav.setAmp(1.)
					accRF_Cav = accRF_Cavs[len(accRF_Cavs) - 1] 
					accRF_Cav.addRF_GapNode(accNode)
					accSeq.addNode(accNode)
				else:
					if(node.getParam("length") != 0.):
						msg = "The LinacLatticeFactory method getLinacAccLattice(names): there is a strange element!"
						msg = msg + os.linesep
						msg = msg + "name=" + node.getName()
						msg = msg + os.linesep
						msg = msg + "type="+node.getType()
						msg = msg + os.linesep
						msg = msg + "length(should be 0.)="+str(node.getParam("length"))
						orbitFinalize(msg)						
					thinNodes.append(node)
			#-----now check the integrety quads and rf_gaps should not overlap
			#-----and create drifts
			copyAccNodes = accSeq.getNodes()[:]
			firstNode = copyAccNodes[0]
			lastNode = copyAccNodes[len(copyAccNodes)-1]
			#insert the drift before the first element if its half length less than its position
			if(math.fabs(firstNode.getLength()/2.0 - firstNode.getParam("pos")) > self.zeroDistance):
				if(firstNode.getLength()/2.0 > firstNode.getParam("pos")):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the first node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + firstNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + firstNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(firstNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(firstNode.getParam("pos"))						
					orbitFinalize(msg)
				else:
					drift = Drift(accSeq.getName()+":drift")
					drift.setLength(firstNode.getParam("pos") - firstNode.getLength()/2.0)
					accSeq.addNode(drift, index = 0)	
			#insert the drift after the last element if its half length less + position is less then the sequence length
			if(math.fabs(lastNode.getLength()/2.0 + lastNode.getParam("pos") - accSeq.getLength()) > self.zeroDistance):
				if(lastNode.getLength()/2.0 + lastNode.getParam("pos") > accSeq.getLength()):
					msg = "The LinacLatticeFactory method getLinacAccLattice(names): the last node is too long!"
					msg = msg + os.linesep
					msg = msg + "name=" + firstNode.getName()
					msg = msg + os.linesep
					msg = msg + "type=" + firstNode.getType()
					msg = msg + os.linesep
					msg = msg + "length=" + str(firstNode.getLength())
					msg = msg + os.linesep
					msg = msg + "pos=" + str(firstNode.getParam("pos"))					
					msg = msg + os.linesep
					msg = msg + "sequence name=" + accSeq.getName()				
					msg = msg + os.linesep
					msg = msg + "sequence length=" + str(accSeq.getLength())			
					orbitFinalize(msg)
				else:
					drift = Drift(accSeq.getName()+":drift")
					drift.setLength(accSeq.getLength() - (lastNode.getParam("pos") + lastNode.getLength()/2.0))
					accSeq.addNode(drift)
			#now move on and generate drifts between (i,i+1) nodes from copyAccNodes
											
			
			
					
						
							
						
		return linacAccLattice
		
		
		



#-----------------------------------------------------
#    LINAC ABSTRACT NODES ELEMENTS
#-----------------------------------------------------

class BaseLinacNode(AccNode):
	""" 
	The base abstract class of the linac accelerator elements hierarchy. 
	It cannot be tilted. The direct subclasses of this class will be markers, 
	user defined nodes etc.
	"""
	def __init__(self, name = "none"):
		"""
		Constructor. Creates the base linac element. This is a superclass for all linac elements.
		"""
		AccNode.__init__(self,name)
		self.setType("baseLinacNode")
		self.__linacSeqence = None
		
	def isRFGap(self):
		"""
		Returns False. The RF Gap node returns True.
		"""
		return False
		
	def setSequence(self, seq):
		"""
		Sets the seqence.
		"""
		self.__linacSeqence = seq
		
	def getSequence(self):
		"""
		Returns the seqence.
		"""
		return self.__linacSeqence
		
	def trackBunch(self, bunch, paramsDict = {}, actionContainer = None):
		"""
		It tracks the bunch through the BaseLinacNode instance.
		"""
		if(actionContainer == None): actionContainer = AccActionsContainer("Linac Bunch Tracking")
		paramsDict["bunch"] = bunch
		
		def track(paramsDict):
			node = paramsDict["node"]
			node.track(paramsDict)
			
		actionContainer.addAction(track, AccActionsContainer.BODY)
		self.trackActions(actionContainer,paramsDict)
		actionContainer.removeAction(track, AccActionsContainer.BODY)		
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the element. Each element 
		should implement this method.
		"""
		pass	

	def trackDesignBunch(self, paramsDict):
		"""
		The RF First Gap nodes will reload this method to setup the design time of passage 
		of the bunch through this node
		"""
		self.track(paramsDict)


class LinacNode(BaseLinacNode):
	"""
	The abstract class of the linac accelerator elements hierarchy that can be tilted.
	"""
	def __init__(self, name = "none"):
		BaseLinacNode.__init__(self,name = "none")
		self.__tiltNodeIN  = TiltElement()
		self.__tiltNodeOUT = TiltElement()
		self.__tiltNodeIN.setName(name+"_tilt_in")
		self.__tiltNodeOUT.setName(name+"_tilt_out")
		self.addChildNode(self.__tiltNodeIN,AccNode.ENTRANCE)
		self.addChildNode(self.__tiltNodeOUT,AccNode.EXIT)
		self.addParam("tilt",self.__tiltNodeIN.getTiltAngle())
		self.setType("linacNode")

	def setTiltAngle(self, angle = 0.):
		"""
		Sets the tilt angle for the tilt operation.
		"""
		self.__params["tilt"] = angle
		self.__tiltNodeIN.setTiltAngle(angle)
		self.__tiltNodeOUT.setTiltAngle( (-1.0) * angle)

	def getTiltAngle(self):
		"""
		Returns the tilt angle for the tilt operation.
		"""
		return self.__tiltNodeIN.getTiltAngle()

	def getNodeTiltIN(self):
		"""
		Returns the Tilt Node instance before this node
		"""
		return self.__tiltNodeIN 
 
	def getNodeTiltOUT(self):
		"""
		Returns the  Tilt Node instance after this node
		"""
		return self.__tiltNodeOUT	

class LinacMagnetNode(LinacNode):
	"""
	The abstract class of the linac magnet.
	"""
	def __init__(self, name = "none"):
		LinacNode.__init__(self,name = "none")
		self.__fringeFieldIN = FringeField(self)
		self.__fringeFieldOUT = FringeField(self)
		self.__fringeFieldIN.setName(name+"_fringe_in")
		self.__fringeFieldOUT	.setName(name+"_fringe_out")	
		self.addChildNode(self.__fringeFieldIN,AccNode.ENTRANCE)
		self.getChildNodes(AccNode.EXIT).insert(0,self.__fringeFieldOUT)
		self.setType("linacMagnet")		
		
	def getNodeFringeFieldIN(self):
		"""
		Returns the FringeField instance before this node
		"""
		return self.__fringeFieldIN 
 
	def getNodeFringeFieldOUT(self):
		"""
		Returns the FringeField instance after this node
		"""
		return self.__fringeFieldOUT 
		
	def setFringeFieldFunctionIN(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch
		through the fringe at the entrance of the node.
		"""
		self.__fringeFieldIN.setFringeFieldFunction(trackFunction)

	def setFringeFieldFunctionOUT(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch
		through the fringe at the exit of the element.
		"""
		self.__fringeFieldOUT.setFringeFieldFunction(trackFunction)

	def getFringeFieldFunctionIN(self, trackFunction = None):
		"""
		Returns the fringe field function that will track the bunch
		through the fringe at the entrance of the node.
		"""
		return self.__fringeFieldIN.getFringeFieldFunction()

	def getFringeFieldFunctionOUT(self, trackFunction = None):
		"""
		Returns the fringe field function that will track the bunch
		through the fringe at the exit of the element.
		"""
		return self.__fringeFieldOUT.getFringeFieldFunction()

	def setUsageFringeFieldIN(self,usage = True):
		"""
		Sets the property describing if the IN fringe
		field will be used in calculation.
		"""
		self.__fringeFieldIN.setUsage(usage)

	def getUsageFringeFieldIN(self):
		"""
		Returns the property describing if the IN fringe
		field will be used in calculation.
		"""
		return self.__fringeFieldIN.getUsage()

	def setUsageFringeFieldOUT(self,usage = True):
		"""
		Sets the property describing if the OUT fringe
		field will be used in calculation.
		"""
		self.__fringeFieldOUT.setUsage(usage)

	def getUsageFringeFieldOUT(self):
		"""
		Returns the property describing if the OUT fringe
		field will be used in calculation.
		"""
		return self.__fringeFieldOUT.getUsage()
		

class TiltNode(BaseLinacNode):
	"""
	The class to do tilt at the entrance of an node.
	"""
	def __init__(self, name = "tilt", angle = 0.):
		"""
		Constructor. Creates the Tilt Node.
		"""
		AccNode.__init__(self,name)
		self.__angle = angle
		self.setType("tilt")

	def setTiltAngle(self, angle = 0.):
		"""
		Sets the tilt angle for the tilt operation.
		"""
		self.__angle = angle

	def getTiltAngle(self):
		"""
		Returns the tilt angle for the tilt operation.
		"""
		return self.__angle

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the titlt node.
		"""
		if(self.__angle != 0.):
			bunch = paramsDict["bunch"]
			TPB.rotatexy(bunch,self.__angle)

class MagnetFringeField(BaseLinacNode):
	"""
	The class is a base class for the fringe field classes for others magnet elements.
	"""
	def __init__(self,  parentNode,  trackFunction = None , name = "fringe field"):
		"""
		Constructor. Creates the Magnet Fringe Field Node.
		"""
		AccNode.__init__(self,name)
		self.setParamsDict(parentNode.getParamsDict())
		self.__trackFunc = trackFunction
		self.__usage = True
		self.setType("magnetFringeField")

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the fringe field node.
		"""
		if(self.__trackFunc != None):
			self.__trackFunc(self,paramsDict)

	def setFringeFieldFunction(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch through the fringe.
		"""
		self.__trackFunc = trackFunction

	def getFringeFieldFunction(self):
		"""
		Returns the fringe field function.
		"""
		return self.__trackFunc

	def setUsage(self,usage = True):
		"""
		Sets the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		self.__usage = usage

	def getUsage(self):
		"""
		Returns the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		return self.__usage

#-----------------------------------------------------
#    LINAC ABST NODES ELEMENTS
#-----------------------------------------------------

class Drift(BaseLinacNode):
	"""
	Drift element.
	"""
	def __init__(self, name = "drift"):
		"""
		Constructor. Creates the Drift element.
		"""
		BaseLinacNode.__init__(self,name)
		self.setType("drift")

	def track(self, paramsDict):
		"""
		The drift class implementation of the AccNode class track(probe) method.
		"""
		length = self.getLength(self.getActivePartIndex())
		bunch = paramsDict["bunch"]
		TPB.drift(bunch, length)
		
class Quad(LinacMagnetNode):
	"""
	Quad Combined Function TEAPOT element.
	"""
	def __init__(self, name = "quad"):
		"""
		Constructor. Creates the Quad Combined Function element .
		"""
		LinacMagnetNode.__init__(self,name)	
		self.addParam("dB/dr",0.)
		self.addParam("poles",[])
		self.addParam("kls",[])
		self.addParam("skews",[])
		self.setnParts(2)
		self.setType("quad")
		
		# B*rho = 3.335640952*momentum [T*m] if momentum in GeV/c
		def fringeIN(node,paramsDict):
			usageIN = node.getUsage()		
			if(not usageIN):
				return
			bunch = paramsDict["bunch"]	
			momentum = bunch.getSyncParticle().momentum()
			kq = node.getParam("dB/dr")/(3.335640952*momentum)
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			TPB.quadfringeIN(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				k = klArr[i]*kq
				skew = skewArr[i]
				TPB.multpfringeIN(bunch,pole,k,skew)

		def fringeOUT(node,paramsDict):
			usageOUT = node.getUsage()
			if(not usageOUT):
				return
			bunch = paramsDict["bunch"]
			momentum = bunch.getSyncParticle().momentum()
			kq = node.getParam("dB/dr")/(3.335640952*momentum)	
			poleArr = node.getParam("poles")
			klArr = node.getParam("kls")
			skewArr = node.getParam("skews")
			length = paramsDict["parentNode"].getLength()
			bunch = paramsDict["bunch"]
			TPB.quadfringeOUT(bunch,kq)
			if(length == 0.):
				return
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				k = klArr[i]*kq
				skew = skewArr[i]
				TPB.multpfringeOUT(bunch,pole,k,skew)

		self.setFringeFieldFunctionIN(fringeIN)
		self.setFringeFieldFunctionOUT(fringeOUT)
		self.getNodeTiltIN().setType("quad tilt in")
		self.getNodeTiltOUT().setType("quad tilt out")
		self.getNodeFringeFieldIN().setType("quad fringe in")
		self.getNodeFringeFieldOUT().setType("quad fringe out")

	def initialize(self):
		"""
		The  Quad Combined Function class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts < 2 and nParts%2 != 0):
			msg = "The Quad Combined Function class instance should have no less than 2 and even number of parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			orbitFinalize(msg)
		lengthIN = (self.getLength()/(nParts - 1))/2.0
		lengthOUT = (self.getLength()/(nParts - 1))/2.0
		lengthStep = lengthIN + lengthOUT
		self.setLength(lengthIN,0)
		self.setLength(lengthOUT,nParts - 1)
		for i in xrange(nParts-2):
			self.setLength(lengthStep,i+1)


	def track(self, paramsDict):
		"""
		The Quad Combined Function TEAPOT  class implementation
		of the AccNode class track(probe) method.
		"""
		bunch = paramsDict["bunch"]		
		momentum = bunch.getSyncParticle().momentum()
		kq = node.getParam("dB/dr")/(3.33564*momentum)		
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getLength(index)
		poleArr = self.getParam("poles")
		klArr = self.getParam("kls")
		skewArr = self.getParam("skews")
		if(index == 0):
			TPB.quad1(bunch, length, kq)
			return
		if(index > 0 and index < (nParts-1)):
			TPB.quad2(bunch, length/2.0)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.quad2(bunch, length/2.0)
			TPB.quad1(bunch, length, kq)
			return
		if(index == (nParts-1)):
			TPB.quad2(bunch, length)
			for i in xrange(len(poleArr)):
				pole = poleArr[i]
				kl = klArr[i]*kq*length/(nParts - 1)
				skew = skewArr[i]
				TPB.multp(bunch,pole,kl,skew)
			TPB.quad2(bunch, length)
			TPB.quad1(bunch, length, kq)
		return		
		
class DCorrectorH(LinacMagnetNode):
	"""
	The Horizontal Dipole Corrector.
	"""
	def __init__(self, name = "correctorh"):
		"""
		Constructor. Creates the Horizontal Dipole Corrector element .
		"""
		LinacMagnetNode.__init__(self,name)
		self.addParam("B",0.)
		self.addParam("effLength",0.)
		self.setType("dch")	
		self.setnParts(1)

	def track(self, paramsDict):
		"""
		The Horizontal Dipole Corrector class implementation of
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getParam("effLength")/nParts
		field = self.getParam("B")
		bunch = paramsDict["bunch"]
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		q = bunch.charge()
		# dp/p = Q*c*B*L/p p in GeV/c c = 2.99792*10^8/10^9		
		kick = q*field*length*0.299792/momentum
		TPB.kick(bunch,kick,0.,0.)

class DCorrectorV(LinacMagnetNode):
	"""
	The Vertical Dipole Corrector.
	"""
	def __init__(self, name = "correctorv"):
		"""
		Constructor. Creates the Vertical Dipole Corrector element .
		"""
		LinacMagnetNode.__init__(self,name)
		self.addParam("B",0.)
		self.addParam("effLength",0.)
		self.setType("dcv")	
		self.setnParts(1)

	def track(self, paramsDict):
		"""
		The Vertical Dipole Corrector class implementation of
		the AccNode class track(probe) method.
		"""
		nParts = self.getnParts()
		index = self.getActivePartIndex()
		length = self.getParam("effLength")/nParts
		field = self.getParam("B")
		bunch = paramsDict["bunch"]
		syncPart = bunch.getSyncParticle()
		momentum = syncPart.momentum()
		q = bunch.charge()
		# dp/p = Q*c*B*L/p p in GeV/c, c = 2.99792*10^8/10^9		
		kick = q*field*length*0.299792/momentum
		TPB.kick(bunch,0,kick,0.)

class BaseRF_Gap(BaseLinacNode):
	"""
	The simplest RF gap representation. The only E*T*L defines all effects of the node.
	"""
	def __init__(self, name = "baserfgap"):
		"""
		Constructor for the simplest RF gap. E0TL parameter is in GeV. Phases are in radians.
		It has 3 parts with lengthes: 0.5 + 0. + 0.5 
		"""
		BaseLinacNode.__init__(self,name)
		self.addParam("E0TL",0.)
		self.addParam("rfCavity", None)
		self.setType("baserfgap")	
		self.__isFirstGap = False
		self.setnParts(3)	

	def initialize(self):
		"""
		The Ring RF TEAPOT class implementation
		of the AccNode class initialize() method.
		"""
		nParts = self.getnParts()
		if(nParts != 3):
			msg = "The simple Rf gap should have 3 parts!"
			msg = msg + os.linesep
			msg = msg + "Method initialize():"
			msg = msg + os.linesep
			msg = msg + "Name of element=" + self.getName()
			msg = msg + os.linesep
			msg = msg + "Type of element=" + self.getType()
			msg = msg + os.linesep
			msg = msg + "nParts =" + str(nParts)
			msg = msg + os.linesep
			msg = msg + "lenght =" + str(self.getLength())
			orbitFinalize(msg)
		length = self.getLength()
		self.setLength(length/2.0,0)
		self.setLength(0.,1)
		self.setLength(length/2.0,2)
	
	def isRFGap(self):
		"""
		Returns True.
		"""
		return True

	def isFirstRFGap(self):
		"""
		Returns True if it is the first gap in RF cavity. 
		"""
		return self.__isFirstGap

	def setAsFirstRFGap(self, isFirst):
		"""
		Sets if it is the first gap in RF cavity. 
		"""
		self.__isFirstGap = isFirst
	
	def setRF_Cavity(self, rf_cav):
		"""
		Sets the parent RF Cavity.
		"""
		self.addParam("rfCavity",rf_cav)

	def getRF_Cavity(self, rf_cav):
		"""
		Returns the parent RF Cavity.
		"""
		return self.getParam("rfCavity")
	
	def track(self, paramsDict):
		"""
		The simplest RF gap class implementation of
		the AccNode class track(probe) method.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		bunch = paramsDict["bunch"]
		if(index == 0 or index == 2):
			TPB.drift(bunch, length)
			return
		E0TL = self.getParam("E0TL")			
		rfCavity = self.getParam("rfCavity")
		frequency = rfCavity.frequency()	
		rfPhase = rfCavity.getPhase()
		phase = rfPhase
		syncPart = bunch.getSyncParticle()
		if(self.__isFirstGap):
			arrival_time = bunch.getSyncParticle().time()
			rfCavity.setFirstGapTime(arrival_time)
			if(rfCavity.isDesignSetUp()):
				designArrivalTime = rfCavity.getDesignArrivalTime()
				phase = math.fmod(frequency*(arrival_time - designArrivalTime)*2.0*math.pi + rfPhase,2.0*math.pi)
		else:
			first_gap_arr_time = rfCavity.getFirstGapTime()
			arrival_time = syncPart.time()
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+rfPhase,2.0*math.pi)	
		#------------------------------------------------------
		# ???? call rf gap with E0TL phase phase of the gap and a longitudinal shift parameter	
		eKin = syncPart.kinEnergy()
		eKin = eKin + E0TL*math.cos(phase)
		syncPart.kinEnergy(eKin)
		
	def trackDesignBunch(self, paramsDict):
		"""
		The RF First Gap node setups the design time of passage 
		of the bunch through this node.
		"""
		index = self.getActivePartIndex()
		length = self.getLength(index)
		bunch = paramsDict["bunch"]
		if(index == 0 or index == 2):
			TPB.drift(bunch, length)
			return		
		E0TL = self.getParam("E0TL")			
		rfCavity = self.getParam("rfCavity")
		frequency = rfCavity.frequency()	
		rfPhase = rfCavity.getDesignPhase()
		phase = rfPhase
		if(self.__isFirstGap):
			arrival_time = bunch.getSyncParticle().time()
			rfCavity.setDesignArrivalTime(arrival_time)
		else:
			first_gap_arr_time = rfCavity.getFirstGapTime()
			arrival_time = bunch.getSyncParticle().time()
			phase = math.fmod(frequency*(arrival_time - first_gap_arr_time)*2.0*math.pi+rfPhase,2.0*math.pi)				
		#------------------------------------------------------
		# ???? call rf gap with E0TL phase phase of the gap and a longitudinal shift parameter	
		eKin = syncPart.kinEnergy()
		eKin = eKin + E0TL*math.cos(phase)
		syncPart.kinEnergy(eKin)	


class TiltElement(BaseLinacNode):
	"""
	The class to do tilt at the entrance of an element.
	"""
	def __init__(self, name = "tilt no name", angle = 0.):
		"""
		Constructor. Creates the Tilt element.
		"""
		AccNode.__init__(self,name)
		self.__angle = angle
		self.setType("tilt teapot")

	def setTiltAngle(self, angle = 0.):
		"""
		Sets the tilt angle for the tilt operation.
		"""
		self.__angle = angle

	def getTiltAngle(self):
		"""
		Returns the tilt angle for the tilt operation.
		"""
		return self.__angle

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the titlt node.
		"""
		if(self.__angle != 0.):
			bunch = paramsDict["bunch"]
			TPB.rotatexy(bunch,self.__angle)

class FringeField(BaseLinacNode):
	"""
	The class is a base class for the fringe field classes for others elements.
	"""
	def __init__(self,  parentNode,  trackFunction = None , name = "fringe field no name"):
		"""
		Constructor. Creates the Fringe Field element.
		"""
		AccNode.__init__(self,name)
		self.setParamsDict(parentNode.getParamsDict())
		self.__trackFunc = trackFunction
		self.__usage = True
		self.setType("fringeField teapot")

	def track(self, paramsDict):
		"""
		It is tracking the dictionary with parameters through
		the fringe field node.
		"""
		if(self.__trackFunc != None):
			self.__trackFunc(self,paramsDict)

	def setFringeFieldFunction(self, trackFunction = None):
		"""
		Sets the fringe field function that will track the bunch through the fringe.
		"""
		self.__trackFunc = trackFunction

	def getFringeFieldFunction(self):
		"""
		Returns the fringe field function.
		"""
		return self.__trackFunc

	def setUsage(self,usage = True):
		"""
		Sets the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		self.__usage = usage

	def getUsage(self):
		"""
		Returns the boolean flag describing if the fringe
		field will be used in calculation.
		"""
		return self.__usage

#----------------------------------------------------------------
# Classes that are specific for the linac model
#----------------------------------------------------------------

class RF_Cavity(NamedObject,ParamsDictObject):
	"""
	This is the class to keep refernces to the RF Gaps which are BaseLinacNode
	subclasses. This class does not belong to the AccNodes.
	"""
	def __init__(self, name = "none"):
		NamedObject.__init__(self, name)
		ParamsDictObject.__init__(self)
		self.__rfGaps = []
		self.addParam("frequency",0.)
		self.addParam("phase",0.)
		self.addParam("amp",0.)		
		self.addParam("designPhase",0.)
		self.addParam("designAmp",0.)		
		self.addParam("firstGapTime",0.)
		self.addParam("designArrivalTime",0.)
		self.addParam("isDesignSetUp",False)
		
	def setDesignSetUp(self,designOnOf):
		""" Sets the design set up information (yes,no). """
		self.setParam("isDesignSetUp",designOnOf)	

	def istDesignSetUp(self):
		""" Returns the design set up information (yes,no). """
		return self.getParam("isDesignSetUp")	
		
	def setFirstGapTime(self,time):
		""" Sets the arrival time for the first RF gap. """
		self.setParam("firstGapTime",time)
		
	def getFirstGapTime(self,time):
		""" Returns the arrival time for the first RF gap. """
		return self.getParam("firstGapTime")

	def setDesignArrivalTime(self,time):
		""" Sets the design arrival time for the first RF gap. """
		self.setParam("designArrivalTime",time)
		
	def getDesignArrivalTime(self,time):
		""" Returns the design arrival time for the first RF gap. """
		return self.getParam("designArrivalTime")
		
	def setDesignPhase(self,phase):
		""" Sets the design phase for the first RF gap. """
		self.setParam("designPhase",phase)
		
	def getDesignPhase(self,time):
		""" Returns the design phase for the first RF gap. """
		return self.getParam("designPhase")

	def setDesignAmp(self,Amp):
		""" Sets the design Amp for the RF cavity. """
		self.setParam("designAmp",Amp)
		
	def getDesignAmp(self,time):
		""" Returns the design Amp for the RF cavity. """
		return self.getParam("designAmp")

	def setPhase(self,phase):
		""" Sets the phase for the first RF gap. """
		self.setParam("phase",phase)
		
	def getPhase(self,time):
		""" Returns the phase for the first RF gap. """
		return self.getParam("phase")

	def setAmp(self,Amp):
		""" Sets the Amp for RF cavity. """
		self.setParam("Amp",Amp)
		
	def getAmp(self,time):
		""" Returns the Amp for RF cavity. """
		return self.getParam("Amp")
		
	def setFrequency(self,freq):
		""" Sets the frequency in Hz. """
		self.setParam("frequency",freq)
		
	def getFrequency(self,time):
		""" Returns the frequency in Hz. """
		return self.getParam("frequency")
		
	def addRF_GapNode(self,rfGap):
		""" Adds the rf gap to the cavity."""
		self.__rfGaps.append(rfGap)
		rfGap.setRF_Cavity(self)
		if(len(self.__rfGaps) == 1):
			rfGap.setAsFirstRFGap(True)
		else:
			rfGap.setAsFirstRFGap(False)
		
	def getRF_GapNodes(self):
		""" Returns the array with rf gaps. """
		return self.__rfGaps[:]
	
class Sequence(NamedObject,ParamsDictObject):
	"""
	This is the class to keep refernces to AccNodes that constitute the accelerator sequence.
	"""
	def __init__(self, name = "none"):
		NamedObject.__init__(self, name)
		ParamsDictObject.__init__(self)
		self.__linacNodes = []
		self.addParam("position",0.)	
		self.addParam("length",0.)	
		
	def addNode(self,node, index = -1):
		""" Adds the Linac Node to the sequence. """
		node.setSequence(self)		
		if(index < 0):
			self.__linacNodes.append(node)
		else:
			self.__linacNodes.insert(index,node)
		
	def getNodes(self):
		""" Returns the array with Linac Nodes. """
		return self.__linacNodes
		
	def setPosition(self, pos):
		""" Sets the position of the sequence. """
		return self.setParam("position",pos)		
		
	def getPosition(self):
		""" Returns the position of the sequence. """
		return self.getParam("position")

	def getLength(self):
		"""
		Returns the total length of the sequence [m].
		"""		
		return self.getParam("length")
		
	def setLength(self, length):
		"""
		Sets the total length of the sequence [m].
		"""		
		return self.setParam("length",length)	
