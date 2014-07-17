#!/usr/bin/env python

#--------------------------------------------------------
# The classes will read data file with beam sizes and transform
# it to another file with less elements (bigger position step)
#--------------------------------------------------------

import math
import sys
import os


mass = 0.939294*1000.

fl_in = open("envelope.txt","r")
lns = fl_in.readlines()
fl_in.close()

fl_out = open("trace_win_results.dat","w")
st = "pos[m] eKin rmsX rmsXp rmsY rmsYp rmsZ rmsdp_p  rmsZp phase time eKin x xp y yp z dp_p zp phase0 time0 eKin0 " 
fl_out.write(st+"\n")

pos_step = 0.05
pos_old = 0.
count = 0
for i in range(1,len(lns)):
	ln = lns[i]
	res_arr = ln.split()
	if(len(res_arr) < 1): continue 
	pos = float(res_arr[0])
	if(pos > (pos_old + pos_step)):
		res_arr = res_arr[0:len(res_arr)-4]
		res_arr[1] = " %12.6f "%(float(res_arr[1])*mass)
		st = ""
		for res in res_arr:
			st += res+" "
		fl_out.write(st+"\n")
		pos_old = pos
		count += 1
fl_out.close()
print "count new lines=",count
	

