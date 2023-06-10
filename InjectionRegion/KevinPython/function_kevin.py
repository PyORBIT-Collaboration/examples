import random
import math

class notRandom:

    def __init__ (self,pos,mom):
        self.name="notRandom"
        self.pos=pos
        self.mom=mom
        
    def getCoordinates(self):
        return (self.pos,self.mom)
        
class xArc:
    def __init__ (self,rmin,rmax,thetaMin,thetaMax,mom):
        self.name="xArc"
        self.rmin=rmin
        self.rmax=rmax
        self.thetaMin=thetaMin
        self.thetaMax=thetaMax
        self.mom=mom
        
    def getCoordinates(self):
    	a=random.random()
    	b=random.random()
    	radius=self.rmin+(self.rmax-self.rmin)*a
    	theta=self.thetaMin+(self.thetaMax-self.thetaMin)*b
    	x=radius*math.cos(theta)
        return (x,self.mom)
	
class yArc:
    def __init__ (self,rmin,rmax,thetaMin,thetaMax,mom):
        self.name="yArc"
        self.rmin=rmin
        self.rmax=rmax
        self.thetaMin=thetaMin
        self.thetaMax=thetaMax
        self.mom=mom
        
    def getCoordinates(self):
    	a=random.random()
    	b=random.random()
    	radius=self.rmin+(self.rmax-self.rmin)*a
    	theta=self.thetaMin+(self.thetaMax-self.thetaMin)*b
    	y=radius*math.sin(theta)
        return (y,self.mom)
        
class zUni:
    def __init__ (self,zmin,zmax,mom):
        self.name="zUni"
        self.zmin=zmin
        self.zmax=zmax
        self.mom=mom        
    def getCoordinates(self):
    	a=random.random()
    	z=self.zmin+(self.zmax-self.zmin)*a
        return (z,self.mom)        