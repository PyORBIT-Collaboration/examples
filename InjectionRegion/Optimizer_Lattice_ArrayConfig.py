#hi

class OptimizerLattice():
	def __init__(self,teapot_lattice):
		self.teapot_lattice=teapot_lattice
		self.doDipoleStrippers=False
		self.dipoleIsStripper=[False,False]
		self.dipoleFixedStripLength=[-1,-1]
		self.dipoleInChicane=[False,False]
		self.dipoleNode=[-1,-1]
		#[0-3] are chicane10-13, 
		self.chicaneNodes=[]
		#[0] is drift DB12 and [1] is drift DB23
		self.driftNodes=[]
		self.chicaneNodeStrength=[[0.041456],[-0.052434],[-0.0298523],[0.0398609]]

	def getDoDipoleStrippers(self):
		return self.doDipoleStrippers
	def setDoDipoleStrippers(self,value):
		self.doDipoleStrippers=value	
	def getDipoleInChicane(self,index):
		return self.dipoleInChicane[index]
	def setDipoleInChicane(self,index,value):
		self.dipoleInChicane[index]=value	
		
	def getDipoleNode(self,index):
		return self.dipoleNode[index]
	def setDipoleNode(self,index,value):
		self.dipoleNode[index]=value	
	
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
		
	def getDipoleIsStripper(self,index):
		return self.dipoleIsStripper[index]
	def setFirstDipoleIsStripper(self,index,value):
		self.dipoleIsStripper[index]=value	
	
	def getDipoleFixedStripLength(self,index):
		return self.dipoleFixedStripLength[index]
	def setFirstDipoleFixedStripLength(self,index,value):
		self.dipoleFixedStripLength[index]=value		
		