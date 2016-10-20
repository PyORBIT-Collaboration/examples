In this directory we have XML files describing the structure of the SNS linac.
The files went through several transformations. The latest one was performed 
on 2016.10.20. We removed all "effLength" parameters from the QUADs and BENDs, 
and specified "length" equal to that value of "effLength". Earlier all quads'
lengths in the linac lattices were defined by the "effLength" inside the linac 
factory. Now we modified the factory to use the "length" parameter.

We also specified the apertures for RF Gaps and BENDS.

We also changed specification of the arrays in XML files from "val," to "val".
The necessary changes were done in the XmlDataAdaptor.py package.

The new files are the following:
sns_linac.xml         - design sns linac with SCL RF up to SCL:23d 
sns_sts_linac.xml     - design secon-target-station sns linac 

