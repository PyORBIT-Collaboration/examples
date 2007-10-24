import sys
import math
import posix

from orbit.lattice import AccLattice, AccElement, AccActionsContainer

lattice = AccLattice("test_lattice")

elem1 = AccElement("el-1")
elem2 = AccElement("el-2")
elem3 = AccElement("el-3")
elem4 = AccElement("el-4")
elem5 = AccElement("el-5")

elem1.setLength(1.1)
elem2.setLength(2.1)
elem3.setLength(3.1)
elem4.setLength(4.1)

lattice.appendChildNode(elem1)
lattice.appendChildNode(elem2)
lattice.appendChildNode(elem3)
lattice.appendChildNode(elem4)
lattice.appendChildNode(elem5)

elem1_1 = AccElement("el-1-1")
elem1_1_1 = AccElement("el-1-1-1")
elem1_1_2 = AccElement("el-1-1-2")
elem1_1_3 = AccElement("el-1-1-3")
elem1_1_4 = AccElement("el-1-1-4")

elem1.appendBodyChildNode(elem1_1)
elem1_1.appendEntranceChildNode(elem1_1_1)
elem1_1.appendExitChildNode(elem1_1_4)
elem1_1.appendBodyChildNode(elem1_1_2)
elem1_1.appendBodyChildNode(elem1_1_3)


elem1_2 = AccElement("el-1-2")
elem2.appendBodyChildNode(elem1_2)

acts = AccActionsContainer()

def Blanks(n):
    s = ""
    for i in xrange(n):
        s += " "
    return s

nLevel = 0

def funcEntrance(paramsDict):
    global nLevel
    nLevel += 1
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"ENTER level=",nLevel," node=",node.getName()

def funcExit(paramsDict):
    global nLevel
    node = paramsDict["node"]
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"EXIT  level=",nLevel," node=",node.getName()
    nLevel -= 1

def funcTrack(paramsDict):
    global nLevel
    node = paramsDict["node"]
    node.track(paramsDict)
    if(paramsDict.has_key("print") and paramsDict["print"] == True):
        print Blanks(nLevel),"TRACK through node =",node.getName()," level=",nLevel

acts.appendEntranceAction(funcEntrance)
acts.appendBodyAction(funcTrack)
acts.appendExitAction(funcExit)

lattice.initialize()

print "Total length=",lattice.getLength()

d = {"print":True}

lattice.trackActions(acts,d)

#========Speed test==========================
count = 0
while(True):
    #lattice.initialize()
    lattice.trackActions(acts)
    if( count % 10000 == 0):
        print "i=",count, " time=",posix.times()[0]
    count += 1

print "====STOP==="

sys.exit(1)


