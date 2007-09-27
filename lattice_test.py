import sys
import math
import posix

from orbit.lattice import AccLattice, AccLine, AccElement, AccActionsContainer

lattice = AccLattice("test_lattice")

line1 = AccLine("line1")
line2 = AccLine("line2")
line3 = AccLine("line3")
line4 = AccLine("line4")

elem1 = AccElement("el-1")
elem2 = AccElement("el-2")
elem3 = AccElement("el-3")
elem4 = AccElement("el-4")

elem5 = AccElement("el-5")

elem1.setLength(1.1)
elem2.setLength(2.1)
elem3.setLength(3.1)
elem4.setLength(4.1)

line1.insertChildNode(elem1)
line2.insertChildNode(elem2)
line3.insertChildNode(elem3)
line4.insertChildNode(elem4)

lattice.insertChildNode(line1)
lattice.insertChildNode(line2)
lattice.insertChildNode(line3)
lattice.insertChildNode(line4)
lattice.insertChildNode(elem5)

elem1_1 = AccElement("el-1-1")
elem1_1_1 = AccElement("el-1-1-1")
elem1_1_2 = AccElement("el-1-1-2")
elem1_1_3 = AccElement("el-1-1-3")
elem1_1_4 = AccElement("el-1-1-4")

elem1.insertChildNodeInBody(elem1_1)
elem1_1.insertChildNodeAtEntrance(elem1_1_1)
elem1_1.insertChildNodeAtExit(elem1_1_4)
elem1_1.insertChildNodeInBody(elem1_1_2)
elem1_1.insertChildNodeInBody(elem1_1_3)


elem1_2 = AccElement("el-1-2")
elem2.insertChildNodeInBody(elem1_2)

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

acts.insertEntranceAction(funcEntrance)
acts.insertBodyAction(funcTrack)
acts.insertExitAction(funcExit)

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


