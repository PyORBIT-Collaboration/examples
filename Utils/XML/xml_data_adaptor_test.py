#-----------------------------------------------------
#XMLDataAdaptor test
#-----------------------------------------------------

import sys

from orbit.utils.xml import XmlDataAdaptor

print "Start."

xml_data_adaptor = XmlDataAdaptor("Test")
xml_data_adaptor.setValue("txt","results")
xml_data_adaptor.setValue("arr_int",[0,1,2,3])
xml_data_adaptor.setValue("arr_double",[-1.,-2.,-3.,-4.])
chld1 = xml_data_adaptor.createChild("child1")
chld2 = chld1.createChild("child2")
chld2.setValue("arr_double1",[-1.,-2.,-3.,-4.])
chld2.setValue("val",0.123456789e-17)
xml_data_adaptor.writeToFile("test.xml")

print "Stop."

