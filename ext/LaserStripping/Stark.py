from starkeffect import *
from ext.las_str.SimplexModMPF import Simplex
from mpmath import *










n1 = 0
n2 = 0
m = 0
point1 = 20


n = n1 + n2 + abs(m) + 1




mp.prec = 100000
Energy = mpf("-0.5")/n/n
Gamma = mpf("0.0")
reZ1 = mpf(2)*(2*n1 + abs(m) + 1)/n
imZ1 = mpf("0.0")
F = mpf("1e-6")


a = Functions(n1, n2, m, point1)
mp.prec = a.setupPrecision(str(F),"("+str(Energy)+" "+ str(Gamma)+")", "("+str(reZ1)+" "+ str(imZ1)+")")


print "mp.prec=  ",mp.dps


#-------------------definition of the M function----------------------------#
def M(args):

   par = a.getM(str(F),"("+str(Energy)+" "+ str(Gamma)+")", "("+str(args[0])+" "+ str(args[1])+")")

   
   line =  par.replace("(","").replace(")","").rsplit(" ")
   comp = mpc(line[0], line[1])
   print "args= ",nstr(args[0],20),nstr(args[1],20),"copm= ",nstr(fabs(comp),20)
   return fabs(comp)
   

#-------------------definition of the M function----------------------------#




#-----------------------Beginning of optimization-----------------------------------------

s = Simplex(M, [reZ1, imZ1], [mpf("1e-5"),mpf("1e-5")])
(values, err, iter) = s.minimize(mpf("0.1"), 10000000,0)
#print values
print err
print iter




