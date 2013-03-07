#-------------------------------------------------------------------------
# This script reads the SuperFish files and generates tables of 
# the TTF functions and their polynomial interpolations. 
# The Super Fish files are the data files with electric and magnetic fields 
# along the RF cavities.
#--------------------------------------------------------------------------


import sys
import math

import orbit_mpi
from orbit_mpi import mpi_comm
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

from spacecharge import Grid2D

from orbit_utils import Function
from orbit_utils import SplineCH
from orbit_utils import GaussLegendreIntegrator
from orbit.utils.fitting import PolynomialFit
from orbit_utils import Polynomial

from orbit.sns_linac.rf_field_readers import SuperFish_3D_RF_FieldReader, RF_AxisFieldAnalysis



# data_in array includes [name_of_file,zSimmetric,beta_min,beta_max, rf_freq]
# zSimmetric = 1 for the cavities symmetrical in z around z=0
# zSimmetric != 1 means there is no symmetry
# beta_min,beta_max - minimum and maximum relativistic beta
# rf_freq - rf frequency in Hz
data_in_arr = []
data_in_arr.append(["scl_medium_beta_rf_cav_field",0,0.4,0.8,805.0e+6])
data_in_arr.append(["scl_high_beta_rf_cav_field",0,0.4,0.9,805.0e+6])
data_in_arr.append(["mebt_1.5cm_field",1,0.04,0.12,402.5e+6])
data_in_arr.append(["mebt_1.8cm_field",1,0.04,0.12,402.5e+6])

dir_name = "data"
n_poly_order = 3

for [name_in,zSimmetric,beta_min,beta_max,rf_freq] in data_in_arr:
	print "============= data file=",name_in,".dat"
	fReader = SuperFish_3D_RF_FieldReader()
	fReader.readFile(dir_name+"/"+name_in+".dat")
	print "beta min, max = %8.4f %8.4f "%(beta_min,beta_max)
	print " z steps =",fReader.zSteps
	print " r steps =",fReader.rSteps
	print " r min max =",fReader.Rmin,"   ",fReader.Rmax
	print " z min max =",fReader.Zmin,"   ",fReader.Zmax	
	spline = fReader.getAxisEz(zSimmetric)
	spline.dump(name_in+"_axis.dat")
	rf_analysis = RF_AxisFieldAnalysis(spline)
	n_table_points = 100
	rf_analysis.makeTransitTimeTables(beta_min,beta_max,n_table_points,rf_freq)
	rf_analysis.makePlynomialFittings(n_poly_order)	
	rf_analysis.dumpTransitTimeTables(dir_name+"/"+name_in+"_t_tp_s_sp.dat")
	rf_analysis.dumpTTFandFitting(name_in+"_ttf_fiitings.dat")	

print "=========================================="
print "Done."




