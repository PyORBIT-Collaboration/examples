#!/usr/bin/env python

#---------------------------------------------------------
#This test will read SAD file and write the MAD file
# with the same structure and parameters
#---------------------------------------------------------

import sys
import math
from sad_parser import SAD_Parser, SAD_LattElement, SAD_LattLine

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
def quadTranslator(elm):
	"""
	MAD - K1=(1/B*pho)*d(B)/d(x)
	SAD - K1 = K1mad*L
	SAD - ROTATE in DEG, but parser already changed it to RAD
	MAD - TILT in RAD
	SAD - MAD have different sign in K1 !!!!
	"""
	L = elm.getParameter("L")
	K1 = elm.getParameter("K1")
	TILT = 0.
	if(elm.hasParameter("ROTATE")):
		TILT = elm.getParameter("ROTATE")
	elm.getParameters().clear()
	elm.getParameters()["L"] = L
	if(L != 0.):
		elm.getParameters()["K1"] = K1/L
	else:
		elm.getParameters()["K1"] = 0.
	if(TILT != 0.):
		elm.getParameters()["TILT"] = TILT

#The DRIFT Translator
def driftTranslator(elm):
	"""
	MAD and SAD - L= length
	"""
	L = elm.getParameter("L")
	elm.getParameters().clear()
	elm.getParameters()["L"] = L

#The APERT Translator
def apertTranslator(elm):
	"""
	MAD - L - length (0 by default)
	MAD - XSIZE - horizontal half aperture
	MAD - YSIZE - vertical half aperture
	SAD - DX1,DX2 and DY1,DY2 - limits
	"""
	DX1 = elm.getParameter("DX1")
	DX2 = elm.getParameter("DX2")
	DY1 = elm.getParameter("DY1")
	DY2 = elm.getParameter("DY2")
	elm.getParameters().clear()
	elm.getParameters()["L"] = 0.
	elm.getParameters()["XSIZE"] = (max((DX1,DX2))-min((DX1,DX2)))/2.0
	elm.getParameters()["YSIZE"] = (max((DY1,DY2))-min((DY1,DY2)))/2.0

#The BEND Translator
def bendTranslator(elm):
	"""
	MAD - E1,E2 see mad user gide
	MAD and SAD - L - length
	SAD - AE1,AE2,E1,E2,FINT,K1
	"""
	L = elm.getParameter("L")
	ANGLE = 0.
	if(elm.hasParameter("ANGLE")):
		ANGLE = elm.getParameter("ANGLE")
	K1 = 0.
	if(elm.hasParameter("K1")):
		K1 = elm.getParameter("K1")
	AE1 = 0.
	if(elm.hasParameter("AE1")):
		AE1 = elm.getParameter("AE1")
	AE2 = 0.
	if(elm.hasParameter("AE2")):
		AE2 = elm.getParameter("AE2")
	E1 = 0.
	if(elm.hasParameter("E1")):
		E1 = elm.getParameter("E1")*ANGLE + AE1
	E2 = 0.
	if(elm.hasParameter("E2")):
		E2 = elm.getParameter("E2")*ANGLE + AE2
	TILT = 0.
	if(elm.hasParameter("ROTATE")):
		TILT = elm.getParameter("ROTATE")
	FINT = 0.
	if(elm.hasParameter("FINT")):
		FINT = elm.getParameter("FINT")
	#put new key-val
	elm.getParameters().clear()
	elm.getParameters()["L"] = L
	if(math.fabs(ANGLE) < 1.0e-6):
		#this will be DRIFT!!!!!
		elm.setType("DRIFT")
		return
	if(K1 != 0. and L != 0.):
		elm.getParameters()["K1"] = K1/L
	if(E1 != 0.):
		elm.getParameters()["E1"] = E1
	if(E2 != 0.):
		elm.getParameters()["E2"] = E2
	elm.getParameters()["ANGLE"] = ANGLE
	if(TILT != 0.):
		elm.getParameters()["TILT"] = TILT
	if(FINT != 0.):
		elm.getParameters()["FINT"] = FINT

#The SEXT Translator
def sextTranslator(elm):
	"""
	MAD - K2,TILT
	MAD and SAD - L - length
	SAD - K2,ROTATE
	"""
	L = elm.getParameter("L")
	K2 = elm.getParameter("K2")
	TILT = 0.
	if(elm.hasParameter("ROTATE")):
		TILT = elm.getParameter("ROTATE")
	#put new key-val
	elm.getParameters().clear()
	elm.getParameters()["L"] = L
	if(L != 0.):
		elm.getParameters()["K2"] = K2/L
	else:
		elm.getParameters()["K2"] = 0.
	if(TILT != 0.):
		elm.getParameters()["TILT"] = TILT

#The MULT Translator
def multTranslator(elm):
	"""
	MAD - up to 9 harmonics, names are KnL and Tn
	MAD and SAD - L - length = 0 by DEFAULT!!!
	SAD - up to 21  harmonics, but higher will be ignored
	SAD - names Kn and SKn
	"""
	knList = []
	for i in xrange(9):
		key = "K"+str(i)
		if(elm.hasParameter(key)):
			knList.append((elm.getParameter(key),i))
	sknList = []
	for i in xrange(9):
		key = "SK"+str(i)
		if(elm.hasParameter(key)):
			knList.append((elm.getParameter(key),i))
	#put new key-val
	elm.getParameters().clear()
	for kn in knList:
		(val,i) = kn
		elm.getParameters()["K"+str(i)+"L"] = val
	for skn in sknList:
		(val,i) = skn
		elm.getParameters()["K"+str(i)+"L"] = val
		elm.getParameters()["T"+str(i)] = math.pi/(2*i+2)

#The CAVI Translator
def caviTranslator(elm):
	"""
	MAD and SAD - L - length
	MAD - VOLT in MV
	SAD - VOLT in V
	SAD - has a lot of params, in MAD all will be ignored
	"""
	L = elm.getParameter("L")
	HARM = elm.getParameter("HARM")
	VOLT = elm.getParameter("VOLT")
	#put new key-val
	elm.getParameters().clear()
	elm.getParameters()["L"] = L
	elm.getParameters()["HARM"] = HARM
	elm.getParameters()["VOLT"] = 0.

#The MARKER Translator
def markerTranslator(elm):
	"""
	MAD and SAD - nothing common
	"""
	elm.getParameters().clear()


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
	translDic = {}
	translDic["QUAD"]  = ("QUADRUPOLE", quadTranslator)
	translDic["DRIFT"] = ("DRIFT", driftTranslator)
	translDic["APERT"] = ("RCOLLIMATOR", apertTranslator)
	translDic["BEND"]  = ("SBEND", bendTranslator)
	translDic["SEXT"]  = ("SEXTUPOLE", sextTranslator)
	translDic["MULT"]  = ("MULTIPOLE", multTranslator)
	translDic["CAVI"]  = ("RFCAVITY", caviTranslator)
	translDic["MARK"]  = ("MARKER", markerTranslator)
	#----------------------------------
	#Translation process
	#----------------------------------
	for elm in elems:
		if(translDic.has_key(elm.getType())):
			(newType,fn) = translDic[elm.getType()]
			elm.setType(newType)
			fn(elm)
		else:
			print "Cannot SAD->MAD translate the type=",elm.getType()
			print "Stop."
			sys.exit(1)


#===========================================
# Main code
#===========================================
if( len(sys.argv) != 3 ):
	print "Usage: >python sad_parser_test.py <name of SAD file to read> <name of Lattice MAD file to write>"
	sys.exit(1)

sad_file_name = sys.argv[1]
mad_file_name = sys.argv[2]

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
lineRING = parser.getSAD_LinesDic()["RING"]
elemsRING = lineRING.getElements()

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
mad_elemsDic = {}
mad_elems = []
for elem in elemsRING:
	if(elem.getType() != "MARKER"):
		if(elem.getType() != "RCOLLIMATOR"):
			mad_elemsDic[elem.getName()] = elem
			mad_elems.append(elem)

print "MAD file will have ",len(mad_elemsDic)," elements."

for elem in mad_elemsDic.values():
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


