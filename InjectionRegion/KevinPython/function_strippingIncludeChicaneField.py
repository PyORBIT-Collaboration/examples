import random
import math
from orbit_utils import Function

class probabilityStrippingWithChicane:

    def __init__ (self,functionx,functiony,n,maxXValue,gamma,beta):
        self.name="probabilityStrippingWithChicane"
        self.magneticFieldx=functionx
        self.magneticFieldy=functiony
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
        
        self.deltayp_rigidity= None  
        self.deltay_rigidity= None 
        self.deltayp_m_rigidity= None  
        self.deltay_m_rigidity= None   
        
        self.InverseFunction= None     
    def tau(self,x):
    	    A1=2.47e-6
    	    A2=4.49e9
    	    c=299792458   
    	    if self.magneticFieldx.getY(x) ==0 and self.magneticFieldy.getY(x) ==0:
    	    	    #tau is infinite. ie it cant be stripped
    	    	    return -1
    	    else : 
    	    	    totalMagneticField=math.sqrt(self.magneticFieldx.getY(x)*self.magneticFieldx.getY(x)+self.magneticFieldy.getY(x)*self.magneticFieldy.getY(x))
    	    	    #print "totalMagneticField= ",totalMagneticField
    	    	    return A1/self.gamma/self.beta/c/totalMagneticField*math.exp(A2/self.gamma/self.beta/c/totalMagneticField)   	    
    def computeFunctions(self):
    	    c=299792458  
    	    stepSize=self.maxXValue/self.n
    	    theSum=0;
    	    theSumDeltaxp=0
    	    theSumDeltax=0
    	    theSumDeltaxp_m=0
    	    theSumDeltax_m=0    
	    
    	    theSumDeltayp=0
    	    theSumDeltay=0
    	    theSumDeltayp_m=0
    	    theSumDeltay_m=0  
	    
    	    self.accumlatedSum= Function()
    	    self.CDF= Function()    
    	    self.deltaxp_rigidity=Function()
    	    self.deltax_rigidity=Function()
    	    self.deltaxp_m_rigidity=Function()
    	    self.deltax_m_rigidity=Function()    
	    
    	    self.deltayp_rigidity=Function()
    	    self.deltay_rigidity=Function()
    	    self.deltayp_m_rigidity=Function()
    	    self.deltay_m_rigidity=Function()    	    
	    
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
    	    	    	    theSumDeltaxp=theSumDeltaxp+self.magneticFieldy.getY(x)*stepSize
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
     	    	    	    	    theSumDeltaxp_m=theSumDeltaxp_m+self.magneticFieldy.getY(x_m)*stepSize
     	    	    	    	    theSumDeltax_m=theSumDeltax_m+(theSumDeltaxp_m+temp_theSumDeltaxp_m)/2.*stepSize
			    self.deltaxp_m_rigidity.add(x,theSumDeltaxp_m)
			    self.deltax_m_rigidity.add(x,theSumDeltax_m)
			    
    	    	    	    temp_theSumDeltayp=theSumDeltayp
    	    	    	    theSumDeltayp=theSumDeltayp+self.magneticFieldx.getY(x)*stepSize
    	    	    	    #this is new one
    	    	    	    theSumDeltay=theSumDeltay+(theSumDeltayp+temp_theSumDeltayp)/2.*stepSize
    	    	    	    
    	    	    	    #this was the og 
    	    	    	    #theSumDeltay=theSumDeltay+theSumDeltayp*stepSize
			    self.deltayp_rigidity.add(x,theSumDeltayp)
			    self.deltay_rigidity.add(x,theSumDeltay)
    	    	    	    
    	    		    theSumDeltayp_m=0
    	                    theSumDeltay_m=0  
     	    	    	    for j in range(i,self.n):
     	    	    	    	    y_m = stepSize*j
     	    	    	    	    temp_theSumDeltayp_m=theSumDeltayp_m
     	    	    	    	    theSumDeltayp_m=theSumDeltayp_m+self.magneticFieldx.getY(y_m)*stepSize
     	    	    	    	    theSumDeltay_m=theSumDeltay_m+(theSumDeltayp_m+temp_theSumDeltayp_m)/2.*stepSize
			    self.deltayp_m_rigidity.add(x,theSumDeltayp_m)
			    self.deltay_m_rigidity.add(x,theSumDeltay_m)			    
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
 	    
    def getdeltayp_rigidity(self):
    	    return self.deltayp_rigidity
    def getdeltay_rigidity(self):
    	    return self.deltay_rigidity    	
    def getdeltayp_m_rigidity(self):
    	    return self.deltayp_m_rigidity
    def getdeltay_m_rigidity(self):
    	    return self.deltay_m_rigidity     	    
    def getInverseFunction(self):
    	    return self.InverseFunction
    def getLength(self,atX):
        return (self.NotNormalizedFunction.getY(atX))
        
