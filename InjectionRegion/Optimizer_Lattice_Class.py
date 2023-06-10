#hi

class OptimizerLattice():
	def __init__(self,teapot_lattice):
		self.teapot_lattice=teapot_lattice
		self.doDipoleStrippers=False
		self.firstDipoleIsStripper=False
		self.secondDipoleIsStripper=False
		self.thirdDipoleIsStripper=False
		self.firstDipoleFixedStripLength=-1
		self.secondDipoleFixedStripLength=-1		
		self.thirdDipoleFixedStripLength=-1
		self.firstDipoleInChicane=False
		self.secondDipoleInChicane=False
		self.thirdDipoleInChicane=False
		self.firstDipoleNode=-1
		self.secondDipoleNode=-1	
		self.thirdDipoleNode=-1
		#[0-3] are chicane10-13, 
		self.chicaneNodes=[]
		#[0] is drift DB12 and [1] is drift DB23
		self.driftNodes=[]
		self.chicaneNodeStrength=[[0.041456],[-0.052434],[-0.0298523],[0.0398609]]

	def getDoDipoleStrippers(self):
		return self.doDipoleStrippers
	def setDoDipoleStrippers(self,value):
		self.doDipoleStrippers=value	
	def getFirstDipoleInChicane(self):
		return self.firstDipoleInChicane
	def setFirstDipoleInChicane(self,value):
		self.firstDipoleInChicane=value	
	def getSecondDipoleInChicane(self):
		return self.secondDipoleInChicane
	def setSecondDipoleInChicane(self,value):
		self.secondDipoleInChicane=value	
	def getThirdDipoleInChicane(self):
		return self.thirdDipoleInChicane
	def setThirdDipoleInChicane(self,value):
		self.thirdDipoleInChicane=value	
		
	def getFirstDipoleNode(self):
		return self.firstDipoleNode
	def setFirstDipoleNode(self,value):
		self.firstDipoleNode=value	
	def getSecondDipoleNode(self):
		return self.secondDipoleNode
	def setSecondDipoleNode(self,value):
		self.secondDipoleNode=value	
	def getThirdDipoleNode(self):
		return self.thirdDipoleNode
	def setThirdDipoleNode(self,value):
		self.thirdDipoleNode=value
		
	def getTeapotLattice(self):
		return self.teapot_lattice
	def setTeapotLattice(self,teapot_lattice):
		self.teapot_lattice=teapot_lattice	
		
	def getChicaneNodes(self):
		return self.chicaneNodes
	def setChicaneNodes(self,value):
		self.chicaneNodes=value	
		
	def getDriftNodes(self):
		return self.driftNodes
	def setDriftNodes(self,value):
		self.driftNodes=value			
		
	def getChicaneNodeStrength(self):
		return self.chicaneNodeStrength
	def setChicaneNodeStrength(self,value):
		self.chicaneNodeStrength=value			
		
	def getFirstDipoleIsStripper(self):
		return self.firstDipoleIsStripper
	def setFirstDipoleIsStripper(self,value):
		self.firstDipoleIsStripper=value	
		
	def getSecondDipoleIsStripper(self):
		return self.secondDipoleIsStripper
	def setSecondDipoleIsStripper(self,value):
		self.secondDipoleIsStripper=value

	def getThirdDipoleIsStripper(self):
		return self.thirdDipoleIsStripper
	def setThirdDipoleIsStripper(self,value):
		self.thirdDipoleIsStripper=value
		
	def getFirstDipoleFixedStripLength(self):
		return self.firstDipoleFixedStripLength
	def setFirstDipoleFixedStripLength(self,value):
		self.firstDipoleFixedStripLength=value		
		
	def getSecondDipoleFixedStripLength(self):
		return self.secondDipoleFixedStripLength
	def setSecondDipoleFixedStripLength(self,value):
		self.secondDipoleFixedStripLength=value			
		
	def getThirdDipoleFixedStripLength(self):
		return self.thirdDipoleFixedStripLength
	def setThirdDipoleFixedStripLength(self,value):
		self.thirdDipoleFixedStripLength=value				