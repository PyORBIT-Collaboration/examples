#hi

class ConfigureFileReader():
	def __init__(self,filename):
		self.filename=filename
		self.theDictionary={}
		self.initialize()		
	def initialize(self):
		openedFile=open("%s"%(self.filename),'r')
		lines=openedFile.readlines()
		for line in lines:
			#skip commentted lines and empty lines
			#print len(line.strip())
			if len(line.strip())>0 and"#"!=line.strip()[0] and "=" in line:
				key_values_pair=line.split("=")
				key=key_values_pair[0].strip()
				values_str=key_values_pair[1]
				values_arr=values_str.strip().split(",")
				self.theDictionary[key]=values_arr
				#print self.theDictionary[key]
				
		openedFile.close()
	def printDictionary(self):
		print self.theDictionary
		
	def getValue(self,key):
		if key not in self.theDictionary:
			print "key=",key," does not exist in filename=",self.filename
		else:
			if len(self.theDictionary[key])==1:
				#print self.theDictionary[key][0]
				return self.theDictionary[key][0]
			else:
				print "key=",key," has more than one value did you mean to use getList?"
	def hasKey(self,key):
		if key in self.theDictionary:
			return True
		else:
			return False
	def getArray(self,key):
		if key not in self.theDictionary:
			print "key=",key," does not exist in filename=",self.filename
		else:
			return self.theDictionary[key]
