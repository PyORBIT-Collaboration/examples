import random
import math
from orbit_utils import Function

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
    def probabilityOfSurvival(self,x):
    	    A1=2.47e-6
    	    A2=4.49e9
    	    c=299792458
    	    if self.magneticField.getY(x) ==0:
    	    	    return 1
    	    else :    	    	   
    	    	    return math.exp(-x*self.magneticField.getY(x)/A1*math.exp(-A2/(self.gamma*self.beta*c*self.magneticField.getY(x))))    
    def tau(self,x):
    	    A1=2.47e-6
    	    A2=4.49e9
    	    c=299792458   
    	    if self.magneticField.getY(x) ==0:
    	    	    #tau is infinite. ie it cant be stripped
    	    	    return -1
    	    else :    	    	   
    	    	    return A1/gamma/beta/c/self.magneticField.getY(x)*math.exp(A2/gamma/beta/c/self.magneticField.getY(x))   	    
    def computeFunctions(self):
    	    stepSize=self.maxXValue/self.n
    	    self.notNormalizedFunction= Function()
    	    self.NormalizedFunction= Function()    
    	    tempFunction=Function()
    	    self.InverseFunction=Function()
   	    normalizeValue=1-self.probabilityOfSurvival(self.maxXValue)
   	    print "normalizeValue= %d" %normalizeValue
    	    for i in range(self.n):
    	    	    x = stepSize*i
    	    	    y = self.probabilityOfSurvival(x)
    	    	    #ynorm=y/normalizeValue
    	    	    ynorm=(1-y)/normalizeValue
    	    	    #and ynorm <(1-1e-4)
    	    	    if y>1e-8 and y <(1-1e-9):
    	    	    	    self.notNormalizedFunction.add(x,y)
    	    	    	    self.NormalizedFunction.add(x,ynorm)
    	    	    	    #tempFunction.add(x,1-ynorm)
    	    	    	    tempFunction.add(x,ynorm)
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
        
