#!/usr/bin/env python

#---------------------------------------------------------
#This test will read SAD file and write the MAD file
# with the same structure and parameters
#---------------------------------------------------------

import sys
import math
from orbit.parsers.sad_parser import SAD_Parser, SAD_LattElement, SAD_LattLine

#===========================================
# Definitions of auxiliary functions
#===========================================
def printMAD_String(fl,res):
	#write everything in 70 simbols in the line
	# fl - file to write
	n = len(res)
	nLength = 66
	i_start = 0
	while(i_start < n):
		i_stop = i_start + nLength
		if(i_stop >= n):
			i_stop = n
			fl.write(res[i_start:i_stop] + " \n")
			break
		else:
			i_stop = res[0:i_stop].rfind(" ")
			fl.write(res[i_start:i_stop]+" &\n")
			i_start = i_stop

#The QUAD Translator
def quadTranslator(elem):
	"""
	MAD - K1=(1/B*pho)*d(B)/d(x)
	SAD - K1 = K1mad*L
	SAD - ROTATE in DEG, but parser already changed it to RAD
	MAD - TILT in RAD
	SAD - MAD have different sign in K1 !!!!
	"""
	L = elem.getParameter("L")
	K1 = elem.getParameter("K1")
	TILT = 0.
	if(elem.hasParameter("ROTATE")):
		TILT = elem.getParameter("ROTATE")
	elem.getParameters().clear()
	elem.getParameters()["L"] = L
	if(L != 0.):
		elem.getParameters()["K1"] = K1/L
	else:
		elem.getParameters()["K1"] = 0.
	if(TILT != 0.):
		elem.getParameters()["TILT"] = TILT

#The DRIFT Translator
def driftTranslator(elem):
	"""
	MAD and SAD - L= length
	"""
	L = elem.getParameter("L")
	elem.getParameters().clear()
	elem.getParameters()["L"] = L

#The APERT Translator
def apertTranslator(elem):
	"""
	MAD - L - length (0 by default)
	MAD - XSIZE - horizontal half aperture
	MAD - YSIZE - vertical half aperture
	SAD - DX1,DX2 and DY1,DY2 - limits
	"""
	DX1 = elem.getParameter("DX1")
	DX2 = elem.getParameter("DX2")
	DY1 = elem.getParameter("DY1")
	DY2 = elem.getParameter("DY2")
	elem.getParameters().clear()
	elem.getParameters()["L"] = 0.
	elem.getParameters()["XSIZE"] = (max((DX1,DX2))-min((DX1,DX2)))/2.0
	elem.getParameters()["YSIZE"] = (max((DY1,DY2))-min((DY1,DY2)))/2.0

#The BEND Translator
def bendTranslator(elem):
	"""
	MAD - E1,E2 see mad user gide
	MAD and SAD - L - length
	SAD - AE1,AE2,E1,E2,FINT,K1
	"""
	L = elem.getParameter("L")
	ANGLE = 0.
	if(elem.hasParameter("ANGLE")):
		ANGLE = elem.getParameter("ANGLE")
	K1 = 0.
	if(elem.hasParameter("K1")):
		K1 = elem.getParameter("K1")
	AE1 = 0.
	if(elem.hasParameter("AE1")):
		AE1 = elem.getParameter("AE1")
	AE2 = 0.
	if(elem.hasParameter("AE2")):
		AE2 = elem.getParameter("AE2")
	E1 = 0.
	if(elem.hasParameter("E1")):
		E1 = elem.getParameter("E1")*ANGLE + AE1
	E2 = 0.
	if(elem.hasParameter("E2")):
		E2 = elem.getParameter("E2")*ANGLE + AE2
	TILT = 0.
	if(elem.hasParameter("ROTATE")):
		TILT = elem.getParameter("ROTATE")
	FINT = 0.
	if(elem.hasParameter("FINT")):
		FINT = elem.getParameter("FINT")
	#put new key-val
	elem.getParameters().clear()
	elem.getParameters()["L"] = L
	if(math.fabs(ANGLE) < 1.0e-6):
		#this will be DRIFT!!!!!
		elem.setType("DRIFT")
		return
	if(K1 != 0. and L != 0.):
		elem.getParameters()["K1"] = K1/L
	if(E1 != 0.):
		elem.getParameters()["E1"] = E1
	if(E2 != 0.):
		elem.getParameters()["E2"] = E2
	elem.getParameters()["ANGLE"] = ANGLE
	if(TILT != 0.):
		elem.getParameters()["TILT"] = TILT
	if(FINT != 0.):
		elem.getParameters()["FINT"] = FINT

#The SEXT Translator
def sextTranslator(elem):
	"""
	MAD - K2,TILT
	MAD and SAD - L - length
	SAD - K2,ROTATE
	"""
	L = elem.getParameter("L")
	K2 = elem.getParameter("K2")
	TILT = 0.
	if(elem.hasParameter("ROTATE")):
		TILT = elem.getParameter("ROTATE")
	#put new key-val
	elem.getParameters().clear()
	elem.getParameters()["L"] = L
	if(L != 0.):
		elem.getParameters()["K2"] = K2/L
	else:
		elem.getParameters()["K2"] = 0.
	if(TILT != 0.):
		elem.getParameters()["TILT"] = TILT

#The MULT Translator
def multTranslator(elem):
	"""
	MAD - up to 9 harmonics, names are KnL and Tn
	MAD and SAD - L - length = 0 by DEFAULT!!!
	SAD - up to 21  harmonics, but higher will be ignored
	SAD - names Kn and SKn
	If L > 0 in SAD it will be a drift in MAD
	"""
	L = 0.
	if(elem.hasParameter("L")):
		L = elem.getParameters()["L"]
	if( L > 0.):
		elem.setType("DRIFT")
		elem.getParameters().clear()
	knList = []
	for i in xrange(9):
		key = "K"+str(i)
		if(elem.hasParameter(key)):
			knList.append((elem.getParameter(key),i))
	sknList = []
	for i in xrange(9):
		key = "SK"+str(i)
		if(elem.hasParameter(key)):
			knList.append((elem.getParameter(key),i))
	#put new key-val
	elem.getParameters().clear()
	if(L > 0.):
		elem.getParameters()["L"] = L	
	for kn in knList:
		(val,i) = kn
		elem.getParameters()["K"+str(i)+"L"] = val
	for skn in sknList:
		(val,i) = skn
		elem.getParameters()["K"+str(i)+"L"] = val
		elem.getParameters()["T"+str(i)] = math.pi/(2*i+2)

#The CAVI Translator
def caviTranslator(elem):
	"""
	MAD and SAD - L - length
	MAD - VOLT in MV
	SAD - VOLT in V
	SAD - has a lot of params, in MAD all will be ignored
	"""
	L = elem.getParameter("L")
	HARM = elem.getParameter("HARM")
	VOLT = elem.getParameter("VOLT")
	#put new key-val
	elem.getParameters().clear()
	elem.getParameters()["L"] = L
	elem.getParameters()["HARM"] = HARM
	elem.getParameters()["VOLT"] = 0.

#The MARKER Translator
def markerTranslator(elem):
	"""
	MAD and SAD - nothing common
	"""
	elem.getParameters().clear()


def SAD_to_MAD_ElementTranslator(elems):
	"""
	Transform all SAD elements from SAD
	notation to MAD. After that the types
	and (key,value) pairs will change.
	"""
	#----------------------------------
	#translation dictionary for types
	#and translation functions
	#----------------------------------
	translDict = {}
	translDict["QUAD"]  = ("QUADRUPOLE", quadTranslator)
	translDict["DRIFT"] = ("DRIFT", driftTranslator)
	translDict["APERT"] = ("RCOLLIMATOR", apertTranslator)
	translDict["BEND"]  = ("SBEND", bendTranslator)
	translDict["SEXT"]  = ("SEXTUPOLE", sextTranslator)
	translDict["MULT"]  = ("MULTIPOLE", multTranslator)
	translDict["CAVI"]  = ("RFCAVITY", caviTranslator)
	translDict["MARK"]  = ("MARKER", markerTranslator)
	translDict["MONI"]  = ("MARKER", markerTranslator)
	#----------------------------------
	#Translation process
	#----------------------------------
	for elem in elems:
		if(translDict.has_key(elem.getType())):
			(newType,fn) = translDict[elem.getType()]
			elem.setType(newType)
			fn(elem)
		else:
			print "Cannot SAD->MAD translate the type=",elem.getType()
			print "Stop."
			sys.exit(1)


#===========================================
# Main code
#===========================================
if( len(sys.argv) != 3 and len(sys.argv) != 4 ):
	print "Usage: >python ",sys.argv[0]," <name of SAD file to read> <name of MAD file to write> [<name of the line>]"
	print "By default the name of the line will be RING"
	sys.exit(1)

sad_file_name = sys.argv[1]
mad_file_name = sys.argv[2]

line_name = "RING"
if(len(sys.argv) == 4):
	line_name = sys.argv[3]

parser = SAD_Parser()
parser.parse(sad_file_name)

lines = parser.getSAD_Lines()
elems = parser.getSAD_Elements()
variables = parser.getSAD_Variables()

print "================================================"
print "The SAD file includes:"
print "Number of lattice accelerator lines     =",len(lines)
print "Number of lattice accelerator elements  =",len(elems)
print "Number of lattice accelerator variables =",len(variables)
print "================================================"

#===================================================
#Translate elements from SAD to MAD
#===================================================
SAD_to_MAD_ElementTranslator(elems)

mad_file = open(mad_file_name,"w")
#===================================================
#write variables needed for tune fitting in MAD
# They are the K1 parameters in QUADS for 7 families
# QDL,QDN,QDX,  QFL,QFM,QFN,QFX
# k1_values{abs(k1*1000000), [quadElements]}
#===================================================
if(not parser.getSAD_LinesDict().has_key(line_name)):
	print "There is no line:",line_name," in the SAD file:",sad_file_name
	print "Stop."
	sys.exit(1)
lineRING = parser.getSAD_LinesDict()[line_name]

elemsRING = lineRING.getElements()

L = 0.
for elem in elemsRING:
	if(elem.hasParameter("L")):
		L = L + elem.getParameter("L")
print "Length of the ",line_name," line=",L

k1_values = {}
for elem in elemsRING:
	if(elem.getType() == "QUADRUPOLE" and elem.hasParameter("K1")):
		if(elem.getParameter("K1") != 0.):
			key = math.floor(math.fabs(elem.getParameter("K1"))*1000000.)
			if(k1_values.has_key(key)):
				k1_values[key].append(elem)
			else:
				k1_values[key] = []
				k1_values[key].append(elem)

print "=====quads families: ======="
i = 0
for key in k1_values.keys():
	i = i + 1
	print "i=",i," K1=",math.fabs(k1_values[key][0].getParameter("K1")),
	print " nQuads=",len(k1_values[key])," nm=",k1_values[key][0].getName()
	#for quad in k1_values[key]:
	#	print "  name=",quad.getName()

#now we replace numerical value of the parameter by variable names
#and write the name := numerical value in the MAD file
for key in k1_values.keys():
	quads = k1_values[key]
	varName = quads[0].getName()[0:3]+"K1"
	val = math.fabs(quads[0].getParameter("K1"))
	mad_file.write(varName+" := "+str(val)+" \n")
	for quad in quads:
		k1 = quad.getParameter("K1")
		if(k1 >= 0.):
			quad.getParameters()["K1"] = varName
		else:
			quad.getParameters()["K1"] = "-"+varName
#---------------------------------------------------------
#Here it will be shortcut. We will dump only elements that
#belong to RING lattice line and we will remove MARKERs
#----------------------------------------------------------
mad_elemsDict = {}
mad_elems = []
for elem in elemsRING:
	if(elem.getType() != "MARKER"):
		if(elem.getType() != "RCOLLIMATOR"):
			mad_elemsDict[elem.getName()] = elem
			mad_elems.append(elem)

print "MAD file will have ",len(mad_elemsDict)," elements."

for elem in mad_elemsDict.values():
	res = ""
	res = res + elem.getName() + " : " + elem.getType()+" "
	for key,val in elem.getParameters().iteritems():
		res = res + ", " + key + " = " + str(val)+ " "
	printMAD_String(mad_file,res)
#--------------------------------------------------------
# Now we create lines that includes not more then nMaxElems
# and RING will consists of these lines
#---------------------------------------------------------
nMaxElems = 75
n = 1 + len(mad_elems)/nMaxElems
tmp_lines =[]
for i in xrange(n):
	tmp_elems = mad_elems[(i*nMaxElems):((i+1)*nMaxElems)]
	if(len(tmp_elems) > 0):
		line_name = "LINETMP0"+str(i)
		res = line_name +" : Line = ( "
		j = 0
		for elem in tmp_elems:
			j = j + 1
			if(j != 1): res = res +" , "
			res = res + elem.getName()
		printMAD_String(mad_file,res+" ) ")
		tmp_lines.append(line_name)
#write RING line
res = "RING : Line = ( "
i = 0
for line_name in tmp_lines:
		i = i + 1
		if(i != 1): res = res +" , "
		res = res + line_name
printMAD_String(mad_file,res+" ) ")
print "The MAD lattice line RING has ",len(tmp_lines)," sublines and ",
print len(mad_elems)," elements."
#=======================================
mad_file.close()
print "Done."
sys.exit(0)


#---------------------------------------------------------------
#This is a previous versions of the elements and lines dumping
#to a MAD file. It will dump everything.
#---------------------------------------------------------------
#===================================================
#write all elements with (key,val) parameters
#===================================================
for elem in elems:
	res = ""
	res = res + elem.getName() + " : " + elem.getType()+" "
	for key,val in elem.getParameters().iteritems():
		res = res + ", " + key + " = " + str(val)+ " "
	printMAD_String(mad_file,res)
#===================================================
#write lines (some lines are included with - sign)
#we have to create additional lines with backwards order
#===================================================
negativeSuffix = "neg"
negativeLineList = []
for line in lines:
	res = line.getName()+" : Line = ( "
	items = line.getItems()
	i = 0
	for item in items:
		i = i + 1
		sign = line.getDirection(item)
		if(i != 1): res = res +" , "
		res = res + item.getName()
		if(sign == -1):
			#print "the line:",line.getName()," has backwards order for:",item.getName(),
			#print " and will be discarded! nElem=",len(item.getElements())
			negativeLineList.append(item)
			res = res + negativeSuffix
	printMAD_String(mad_file,res+" ) ")
#===================================================
#printing the reverse lines with new names
for line in negativeLineList:
	res = line.getName()+negativeSuffix+" : Line = ( "
	items = line.getElements()
	items.reverse()
	i = 0
	for item in items:
		i = i + 1
		if(i != 1): res = res +" , "
		res = res + item.getName()
	printMAD_String(mad_file,res + " ) ")
#===================================================


