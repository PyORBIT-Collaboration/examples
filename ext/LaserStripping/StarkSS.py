from ext.las_str.StarkEffectSS  import *
from mpmath import *
import os

orbit_path = os.environ["ORBIT_ROOT"]



(n1,n2,m) = (0,0,1)




b = Stark_calcSS(n1, n2, m, 10)





b.F = mpf("6.0e-3")

b.main()

E = b.Energy
#    G = b.readG(orbit_path+"/ext/laserstripping/transitions/")
G = mpf("1.3e-4")
    

for i in range(-10,11):
    Ei = E + i*G/1
    f = open('TDM_%i%i%i_F=%f.txt'%(b.n1,b.n2,b.m,b.F),'a')
    print >>f,nstr(Ei,30),"\t",nstr(b.TDM(Ei),30)
    f.close()
sys.exit()

   
#    f = open('%i%i%i.txt'%(n1,n2,m),'a')
#    print >>f, nstr(b.F,50).ljust(30),nstr(b.Energy,20).ljust(30),nstr(b.calc,20).ljust(30)
#    f.close()
    
#    f = open('_%i%i%i.txt'%(n1,n2,m),'a')
#    print >>f, b.pointM,"\t",b.pointN,"\t",b.calc,"\t",str(b.Energy),"\t",str(b.Gamma),"\t",str(re(b.Z1)),"\t",str(im(b.Z1)),"\t",nstr(b.F,50)
#    f.close()
    