import math
from orbit.teapot.teapot import NodeTEAPOT

class MaxDisplacement_Node(NodeTEAPOT):
    #updates max displacement of bunch
    def __init__(self,name = "MaxDisplacementNode"):
        NodeTEAPOT.__init__(self,name)


    def track(self, paramsDict):
        """
        The WS_AccNode class implementation of the AccNode class track(probe) method.
        """
        bunch = paramsDict["bunch"]
        if(not paramsDict.has_key("maxXDisplacement")): paramsDict["maxXDisplacement"] = 0.
        maxX=paramsDict["maxXDisplacement"]
	for i in range(bunch.getSize()):   
		if math.fabs(bunch.x(i))>math.fabs(paramsDict["maxXDisplacement"]):
			paramsDict["maxXDisplacement"]=bunch.x(i)


