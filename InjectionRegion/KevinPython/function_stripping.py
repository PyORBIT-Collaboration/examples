import random
import math
from orbit.utils import Function

class probabilityStripping:

    def __init__ (self,function,n,maxXValue,gamma,beta):
        self.name="probabilityStripping"
        self.magneticField=function
        self.maxXValue=maxXValue
        self.n=n
        self.gamma=gamma
        self.beta=beta
        self.notNormalizedFunction= None
        self.NormalizedFunction= None
        self.InverseFunction= None    
    def probabilityOfSurvival (self,x):
    	    A1=2.47e-6
    	    A2=4.49e9
    	    c=299792458
    	    return math.exp(-x*self.magneticField.getY(x)/A1*math.exp(-A2/(self.gamma*self.beta*c*self.magneticField.getY(x))))
    def computeFunctions(self):
    	    stepSize=self.maxXValue/self.n
    	    self.notNormalizedFunction= Function()
    	    self.NormalizedFunction= Function()    
    	    tempFunction=Function()
   	    normalizeValue=1-probabilityOfSurvival(self.maxXValue)
    	    for i in range(self.n):
    	    	    x = stepSize*i
    	    	    y = probabilityOfSurvival(x)
    	    	    ynorm=y/normalizeValue
    	    	    self.notNormalizedFunction.add(x,y)
    	    	    self.NormalizedFunction.add(x,ynorm)
    	    	    tempFunction.add(x,1-ynorm)
    	    wasSuccess=tempFunction.setInverse(self.InverseFunction)
    	    if wasSuccess is 1:
    	    	    print "successfully Inverted"
    	    else:
    	    	    print "failed to Invert"
    def getNotNormalizedFunction(self):
    	    return self.notNormalizedFunction
    def getNormalizedFunction(self):
    	    return self.NormalizedFunction
    def getInverseFunction(self):
    	    return self.InverseFunction
    def getLength(self,atX):
        return (self.NotNormalizedFunction.getY(atX))
        
