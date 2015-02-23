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
xml_text = xml_data_adaptor.makeXmlText()
print xml_text
xml_data_adaptor.writeToFile("test.xml")
print "==================writing to a file is done============================"
print "======================================================================="

xml_data_adaptor_new = XmlDataAdaptor.adaptorForFile("test.xml")
xml_text = xml_data_adaptor_new.makeXmlText()
print "====================After reading the file============================="
print xml_text
print "======================================================================="
child1 = xml_data_adaptor_new.childAdaptors("child1")[0]
child2 = child1.childAdaptors("child2")[0]
print " child2 int arr =",child2.intArrayValue("arr_double1")
print " child2 double arr =",child2.doubleArrayValue("arr_double1")
print " child2 double val =",child2.doubleArrayValue("val")
print " child2 double val =",child2.doubleValue("val")

print "Stop."

