import math
from ConfigureFileClass import ConfigureFileReader

class MagneticField():
	def __init__(self,filename,nodeName,isStripper,refNodeName,position,fixedStrippingLength):
		self.filename=filename
		self.magneticFieldDictionary=ConfigureFileReader(self.filename)
		#self.magneticFieldDictionary.printDictionary()	
		self.nParts=-1
		self.stripperLength=-1
		self.stripperStrengthMax=-1
		self.stripperStrengthMin=-1
		self.cutLength=-1		
		self.fieldDirection=0
		self.isStripper=isStripper
		self.fixedStrippingLength=fixedStrippingLength
		self.nodeName=nodeName
		self.nodePosition=position
		self.nodeIndex=-1
		self.refNodeName=refNodeName
		self.refNodeIndex=-1
		self.isInsideChicane=False
		#chicane numbering starts at 0,1,2,3
		self.chicaneItIsInside=-1
		#0 is uniform
		#1 is linear ramp
		self.magnetFunctionType=1
		self.foilTest=False
		self.initialize()
	
	def getStripperLength(self):
		return self.stripperLength
	def setStripperLength(self,value):
		self.stripperLength=value
		
	def getStripperStrengthMax(self):
		return self.stripperStrengthMax
	def setStripperStrengthMax(self,value):
		self.stripperStrengthMax=value
		
	def getStripperStrengthMin(self):
		return self.stripperStrengthMin
	def setStripperStrengthMin(self,value):
		self.stripperStrengthMin=value
		
	def getCutLength(self):
		return self.cutLength
	def setCutLength(self,value):
		self.cutLength=value
		
	def getStripperLength(self):
		return self.stripperLength
	def setStripperLength(self,value):
		self.stripperLength=value
		
	def getFieldDirection(self):
		return self.fieldDirection
	def setFieldDirection(self,value):
		self.fieldDirection=value
	def getIsStripper(self):
		return self.isStripper
	def setIsStripper(self,value):
		self.isStripper=value	
	def getMagnetFunctionType(self):
		return self.magnetFunctionType
	def setMagnetFunctionType(self,value):
		self.magnetFunctionType=value
		
	def getNodeName(self):
		return self.nodeName
	def setNodeName(self,value):
		self.nodeName=value		

	def getNodePosition(self):
		return self.nodePosition
	def setNodePosition(self,value):
		self.nodePosition=value
		
	def getNodeIndex(self):
		return self.nodeIndex
	def setNodeIndex(self,value):
		self.nodeIndex=value
		
	def getRefNodeName(self):
		return self.refNodeName
	def setRefNodeName(self,value):
		self.refNodeName=value	
		
	def getRefNodeIndex(self):
		return self.refNodeIndex
	def setRefNodeIndex(self,value):
		self.refNodeIndex=value
		
	def getFixedStrippingLength(self):
		return self.fixedStrippingLength
	def setFixedStrippingLength(self,value):
		self.fixedStrippingLength=value
		
	def getIsInsideChicane(self):
		return self.isInsideChicane
	def setIsInsideChicane(self,value):
		self.isInsideChicane=value

	def getChicaneItIsInside(self):
		return self.chicaneItIsInside
	def setChicaneItIsInside(self,value):
		self.chicaneItIsInside=value

	def getNParts(self):
		return self.nParts
	def setNParts(self,value):
		self.nParts=value	
		
	def getFoilTest(self):
		return self.foilTest
	def setFoilTest(self,value):
		self.foilTest=value
		
	def getValueOfField(self,value):
		y="temp"
		if self.magnetFunctionType==1:
			if value<self.cutLength :
				y= (self.stripperStrengthMax-self.stripperStrengthMin)/self.cutLength*value+self.stripperStrengthMin
			elif value>=self.cutLength:
				y= self.stripperStrengthMax
		elif self.magnetFunctionType==0:
			y= self.stripperStrengthMax
		else:
			print "self.magnetFunctionType=",self.magnetFunctionType," is not defined"
		return y
	def initialize(self):
		self.stripperLength=float(self.magneticFieldDictionary.getValue("stripperLength"))
		self.stripperStrengthMax=float(self.magneticFieldDictionary.getValue("stripperStrengthMax"))
		self.stripperStrengthMin=float(self.magneticFieldDictionary.getValue("stripperStrengthMin"))
		self.cutLength=float(self.magneticFieldDictionary.getValue("cutLength"))	
		if self.magneticFieldDictionary.getValue("fieldDirection").lower()=="up":
			self.fieldDirection=math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection").lower()=="down":
			self.fieldDirection=-math.pi/2.
		elif self.magneticFieldDictionary.getValue("fieldDirection").lower()=="left":
			self.fieldDirection=0
		elif self.magneticFieldDictionary.getValue("fieldDirection").lower()=="right":
			self.fieldDirection=math.pi
		else:
			self.fieldDirection=float(self.magneticFieldDictionary.getValue("fieldDirection"))			
	
		#self.isStripper=self.magneticFieldDictionary.getValue("isStripper")
		self.magnetFunctionType=int(self.magneticFieldDictionary.getValue("magnetFunction"))
		self.nParts=int(self.magneticFieldDictionary.getValue("nParts"))
		if self.magneticFieldDictionary.hasKey("foilTest"):
			if self.magneticFieldDictionary.getValue("foilTest") == "True":
				self.foilTest=True