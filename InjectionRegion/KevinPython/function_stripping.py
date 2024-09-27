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
        self.accumlatedSum= None
        self.CDF= None
        self.deltaxp_rigidity= None  
        self.deltax_rigidity= None 
        self.deltaxp_m_rigidity= None  
        self.deltax_m_rigidity= None         
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
    	    	    return A1/self.gamma/self.beta/c/self.magneticField.getY(x)*math.exp(A2/self.gamma/self.beta/c/self.magneticField.getY(x))   	    
    def computeFunctions(self):
    	    c=299792458  
    	    stepSize=self.maxXValue/self.n
    	    theSum=0;
    	    theSumDeltaxp=0
    	    theSumDeltax=0
    	    theSumDeltaxp_m=0
    	    theSumDeltax_m=0    	    
    	    self.accumlatedSum= Function()
    	    self.CDF= Function()    
    	    self.deltaxp_rigidity=Function()
    	    self.deltax_rigidity=Function()
    	    self.deltaxp_m_rigidity=Function()
    	    self.deltax_m_rigidity=Function()    	    
    	    self.InverseFunction=Function()
   	    #normalizeValue=1-self.probabilityOfSurvival(self.maxXValue)
   	    #print "normalizeValue= %d" %normalizeValue
    	    for i in range(self.n):
    	    	    x = stepSize*i
    	    	    y = self.tau(x)
    	    	    if y<0:
    	    	    	    #do nothing because tau is infinite
    	    	    	    pass
    	    	    else:
    	    	    	    theSum=theSum+stepSize/self.gamma/self.beta/c/y
    	    	    	    temp_theSumDeltaxp=theSumDeltaxp
    	    	    	    theSumDeltaxp=theSumDeltaxp+self.magneticField.getY(x)*stepSize
    	    	    	    #this is new one
    	    	    	    theSumDeltax=theSumDeltax+(theSumDeltaxp+temp_theSumDeltaxp)/2.*stepSize
    	    	    	    
    	    	    	    #this was the og 
    	    	    	    #theSumDeltax=theSumDeltax+theSumDeltaxp*stepSize
			    self.deltaxp_rigidity.add(x,theSumDeltaxp)
			    self.deltax_rigidity.add(x,theSumDeltax)
    	    	    	    
    	    		    theSumDeltaxp_m=0
    	                    theSumDeltax_m=0  
     	    	    	    for j in range(i,self.n):
     	    	    	    	    x_m = stepSize*j
     	    	    	    	    temp_theSumDeltaxp_m=theSumDeltaxp_m
     	    	    	    	    theSumDeltaxp_m=theSumDeltaxp_m+self.magneticField.getY(x_m)*stepSize
     	    	    	    	    theSumDeltax_m=theSumDeltax_m+(theSumDeltaxp_m+temp_theSumDeltaxp_m)/2.*stepSize
			    self.deltaxp_m_rigidity.add(x,theSumDeltaxp_m)
			    self.deltax_m_rigidity.add(x,theSumDeltax_m)
			    #self.deltax_rigidity.add(x,theSumDeltax+(self.maxXValue-(x+stepSize))*theSumDeltaxp) 
   	    	    	    self.accumlatedSum.add(x,theSum)
    	    	    	    #3e-7 gives 1-exp(-15)=exp(-3e-7)
    	    	    	    if theSum<15 and theSum>3e-7: 
    	    	    	    	    self.CDF.add(x,1-math.exp(-theSum))

    	    	    #ynorm=y/normalizeValue
    	    	    #ynorm=(1-y)/normalizeValue
    	    	    #and ynorm <(1-1e-4)
    	    	    #if y>1e-8 and y <(1-1e-9):
    	    	    	    #self.notNormalizedFunction.add(x,y)
    	    	    	    #self.NormalizedFunction.add(x,ynorm)
    	    	    	    #tempFunction.add(x,1-ynorm)
    	    	    	    #tempFunction.add(x,ynorm)
    	      
    	    wasSuccess=self.CDF.setInverse(self.InverseFunction)
    	    if wasSuccess is 1:
    	    	    print "successfully Inverted"
    	    else:
    	    	    print "failed to Invert"
    def getaccumlatedSum(self):
    	    return self.accumlatedSum
    def getCDF(self):
    	    return self.CDF
    def getdeltaxp_rigidity(self):
    	    return self.deltaxp_rigidity
    def getdeltax_rigidity(self):
    	    return self.deltax_rigidity    	
    def getdeltaxp_m_rigidity(self):
    	    return self.deltaxp_m_rigidity
    def getdeltax_m_rigidity(self):
    	    return self.deltax_m_rigidity     	    
    def getInverseFunction(self):
    	    return self.InverseFunction
    def getLength(self,atX):
        return (self.NotNormalizedFunction.getY(atX))
        
