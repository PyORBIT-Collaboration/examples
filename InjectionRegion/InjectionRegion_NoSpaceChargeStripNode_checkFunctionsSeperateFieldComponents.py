#############################################################
# This script reads the input MAD file with lattice information,
# creates the TEAPOT lattice, and modifies this lattice by inserting 
# injection nodes
##############################################################

from ROOT import *

import math
import sys

from orbit_utils import Function
from KevinPython.function_stripping import probabilityStripping
from KevinPython.function_strippingIncludeChicaneField import probabilityStrippingWithChicane

e_kin_ini = 1.0 # in [GeV]
mass =  0.93827231 #0.939294    # in [GeV]
gamma = (mass + e_kin_ini)/mass
beta = math.sqrt(gamma*gamma - 1.0)/gamma
print "relat. gamma=",gamma
print "relat.  beta=",beta

theEffLength=0.03*2
#theEffLength=0.01
#fieldStrength=1.3
#fieldStrengthMin=.2
fieldStrength=.07
fieldStrengthMin=.07
cutLength=0.03
fieldDirection=math.pi/2.

A1=2.47e-6
A2=4.49e9

def constantField(x):
	return fieldStrength
def pieceWiseField(x):
	if x<=0.0:
		return 0
	elif x<cutLength :
		return fieldStrength/cutLength*x
	elif x>=cutLength:
		return fieldStrength
	pass
def pieceWiseField2(x):
	if x<cutLength :
		return (fieldStrength-fieldStrengthMin)/cutLength*x+fieldStrengthMin
	elif x>=cutLength:
		return fieldStrength
	pass
magneticFieldx= Function()
magneticFieldy= Function()
n=1000
maxValue=theEffLength
step=maxValue/n

for i in range(n):
	x = step*i;
	#y = constantField(x)
	y = pieceWiseField2(x)
	magneticFieldx.add(x,y*math.cos(fieldDirection))
	magneticFieldy.add(x,y*math.sin(fieldDirection))
	
theStrippingFunctions=probabilityStrippingWithChicane(magneticFieldx,magneticFieldy,n,maxValue,gamma,beta)
theStrippingFunctions.computeFunctions()
accumlatedSum=theStrippingFunctions.getaccumlatedSum()
CDF=theStrippingFunctions.getCDF()
deltaxp_rigidity=theStrippingFunctions.getdeltaxp_rigidity()
deltax_rigidity=theStrippingFunctions.getdeltax_rigidity()
deltaxp_m_rigidity=theStrippingFunctions.getdeltaxp_m_rigidity()
deltax_m_rigidity=theStrippingFunctions.getdeltax_m_rigidity()

deltayp_rigidity=theStrippingFunctions.getdeltayp_rigidity()
deltay_rigidity=theStrippingFunctions.getdeltay_rigidity()
deltayp_m_rigidity=theStrippingFunctions.getdeltayp_m_rigidity()
deltay_m_rigidity=theStrippingFunctions.getdeltay_m_rigidity()	

InverseFunction=theStrippingFunctions.getInverseFunction()
  


theFunctions=[accumlatedSum,CDF,InverseFunction,magneticFieldx,magneticFieldy,deltaxp_rigidity,deltax_rigidity,deltayp_rigidity,deltay_rigidity]
titles=["accumlatedSum","CDF","InverseFunction","magneticFieldx","magneticFieldy","deltaxp_rigidity","deltax_rigidity","deltayp_rigidity","deltay_rigidity"]
notNormalizedGraph= TGraph()
NormalizedGraph= TGraph()
InverseGraph= TGraph()
counter=0
for function in theFunctions:
	currentGraph=TGraph()
	for i in range(function.getSize()):
		currentGraph.SetPoint(currentGraph.GetN(),float(function.x(i)),float(function.y(i)))
	theCanvas=TCanvas("TGraph","TGraph",0,0,500,500)
	currentGraph.Draw("AP")   
	theCanvas.Print("%s.png"%(titles[counter]))
	#theCanvas.SaveAs("temp%d.png"%counter)
	theCanvas.Clear()
	counter=counter+1


print "notNormalizedFunction->getY(maxValue)"
#print accumlatedSum.getY(maxValue)
print CDF.getY(maxValue)
print deltaxp_rigidity.getY(maxValue)
print deltax_rigidity.getY(maxValue)
print magneticFieldy.getY(maxValue)
print magneticFieldx.getY(maxValue)