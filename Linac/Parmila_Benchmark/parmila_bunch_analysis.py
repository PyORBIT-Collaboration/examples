#!/usr/bin/env python

#--------------------------------------------------------
# The script will calculate the Twiss parameters of the 
# Parmila particles' distribution that was dumped and 
# transformed to the text file by readdst.exe
#--------------------------------------------------------

import math
import sys

from orbit.bunch_generators import TwissAnalysis

if(len(sys.argv) != 2):  
	print "Usage: >python ",sys.argv[0]," <name of Parmila bunch text file>"  
	sys.exit(1)
	
parmila_bunch_in = open(sys.argv[1],"r")

twiss_analysis = TwissAnalysis(3)


s = parmila_bunch_in.readline()
while(s.find("x(cm)             xpr(=dx/ds)") < 0):
	s = parmila_bunch_in.readline()

s = parmila_bunch_in.readline()
n_parts = 0
while(len(s) > 20):
	s_arr = s.split()
	if(len(s_arr) != 6):
		print "Structure of the Parmila dump file is wrong! file name =",sys.argv[1]
		sys.exit(1)
	n_parts = n_parts + 1
	x = float(s_arr[0])
	xp = float(s_arr[1])
	y = float(s_arr[2])
	yp = float(s_arr[3])
	phi = float(s_arr[4])
	Ek = float(s_arr[5])
	twiss_analysis.account((x,xp,y,yp,phi*180/math.pi,Ek))
	s = parmila_bunch_in.readline()
	
parmila_bunch_in.close()

print "n particles =",n_parts
mass = 939.3014      # MeV - mass of H-
e_kin = twiss_analysis.getAvgA_AP(2)[1]
gamma = 1.0 + e_kin/mass
beta = math.sqrt(1.0 - 1.0/gamma**2)
bg = gamma*beta

#------------------------------
# print the parameters of the distribution
#------------------------------
print "file =",sys.argv[1]
(alpha_x,beta_x,gamma_x,emitt_x) = twiss_analysis.getTwiss(0)
(alpha_y,beta_y,gamma_y,emitt_y) = twiss_analysis.getTwiss(1)
(alpha_z,beta_z,gamma_z,emitt_z) = twiss_analysis.getTwiss(2)
emitt_x = 1000.*emitt_x*bg
emitt_y = 1000.*emitt_y*bg

print "-------------------bunch's twiss parameters----------------------------------------"
print "X alpha= %12.5g    beta [cm/rad]  =%12.5g    gamma = %12.5g  norm. emitt[cm*mrad]  = %12.5g "%(alpha_x,beta_x,gamma_x,emitt_x)
print "Y alpha= %12.5g    beta [cm/rad]  =%12.5g    gamma = %12.5g  norm. emitt[cm*mrad]  = %12.5g "%(alpha_y,beta_y,gamma_y,emitt_y)
print "Z alpha= %12.5g    beta [deg/MeV] =%12.5g    gamma = %12.5g  emitt[deg*MeV] = %12.5g "%(alpha_z,beta_z,gamma_z,emitt_z)

print "-------------------centroid params--------------------------------------"
print "X x_avg [cm]    =%12.5g    xp_avg [rad]  =%12.5g "%twiss_analysis.getAvgA_AP(0)
print "Y y_avg [cm]    =%12.5g    xp_avg [rad]  =%12.5g "%twiss_analysis.getAvgA_AP(1)
print "Z phi_avg [deg] =%12.5g    Ek [Mev]      =%12.5g "%twiss_analysis.getAvgA_AP(2)

(x_rms,xp_rms) =  twiss_analysis.getRmsA_AP(0)
(y_rms,yp_rms) =  twiss_analysis.getRmsA_AP(1)
(z_rms,zp_rms) =  twiss_analysis.getRmsA_AP(2)
print "-------------------Rms--------------------------------------"
print "X x_rms [cm]    =%12.5g  xp_rms [deg] =%12.5g "%(x_rms,xp_rms*180./math.pi)
print "Y y_rms [cm]    =%12.5g  yp_rms [deg] =%12.5g "%(y_rms,yp_rms*180./math.pi)
print "Z phi_rms [deg] =%12.5g  Ek_rms [MeV] =%12.5g "%(z_rms,zp_rms)

(x_max,xp_max) =  twiss_analysis.getMaxA_AP(0)
(x_min,xp_min) =  twiss_analysis.getMinA_AP(0)
(y_max,yp_max) =  twiss_analysis.getMaxA_AP(1)
(y_min,yp_min) =  twiss_analysis.getMinA_AP(1)
(z_max,zp_max) =  twiss_analysis.getMaxA_AP(2)
(z_min,zp_min) =  twiss_analysis.getMinA_AP(2)
print "-------------------Min Max--------------------------------------"
print "X x_min_max [cm]    =%12.5g %12.5g  xp_min_max [rad] =%12.5g %12.5g "%(x_min,x_max,xp_min,xp_max)
print "Y y_min_max [cm]    =%12.5g %12.5g  yp_min_max [rad] =%12.5g %12.5g "%(y_min,y_max,yp_min,yp_max)
print "Z phi_min_max [deg] =%12.5g %12.5g  Ek_min_max [MeV] =%12.5g %12.5g "%(z_min,z_max,zp_min,zp_max)

