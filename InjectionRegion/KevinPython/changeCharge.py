
from orbit.teapot.teapot import NodeTEAPOT

class Change_Charge_Child(NodeTEAPOT):
    """
    Changes charge of bunch
    """
    def __init__(self,name = "ChangeCharge",charge=1):
        NodeTEAPOT.__init__(self,name)
        self.charge=charge

    def track(self, paramsDict):
        """
        The WS_AccNode class implementation of the AccNode class track(probe) method.
        """
        debug=False
        bunch = paramsDict["bunch"]
        bunch.charge(self.charge)
        #print bunch.charge()
 		