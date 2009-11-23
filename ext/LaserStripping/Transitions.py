from ext.las_str.StarkEffect  import *
from mpmath import *
import os, time
orbit_path = os.environ["ORBIT_ROOT"]

def saferead(addr,j):
    while True:
        time.sleep(1)
        a = ['','','','','','','','','','']
        while True:
            for i in xrange(10):
                time.sleep(0.1)
                f = open(addr,'r')
                a[i] = f.read()
                f.close()

            if(a[0]==a[1]==a[2]==a[3]==a[4]==a[5]==a[6]==a[7]==a[8]==a[9]):
                break
        if(len(a[0].splitlines())>j):
            break
    return a[0].splitlines()[j]




"""
Energy =  "-0.078388135826455327502725226063798950806113275840260409954515146743327933756257367540920760951273929093889694278127628746102"   
Gamma =  "0.0040161803800431632763025263340667207336114464420696927055098670080233581163098292671367803013726018619806251669267026763723"

reZ1 =  "0.81585842928137492019582005646046483363587645751078955361133176141506139098267257422827215926748915548775364335273016790713"
imZ1 =  "0.0095842521743217275590859123877076904119081589055758509138914523361076743965651584238851399334901989887306579877705173570832"

F = "0.002"

maxM = 18
maxN = 46


a = WaveFunc(0,2,0,maxM,maxN,Energy,Gamma,reZ1,imZ1,F)


z0 = 1/sqrt(mpf(F))

print "mode = ", a.mode
M = a.M
N = a.N



mp.dps = 30
print "integration"
print quad(lambda mu: quad(lambda nu: fabs(M(mu)*N(nu))**2,[0,sqrt(mu*mu+2*z0)]),[0,maxM])
"""



#sqrt_norm=  1.77245508719921



(n1,n2,m) = (0,1,1)

num0 = 0
num = 0


for k in xrange(0,1000000000):
#    if not os.path.exists("_run"):    sys.exit(0)

    ground = saferead('_000.txt',k)
    upper = saferead('_%i%i%i.txt'%(n1,n2,m),k)
    
    F = ground.split()[7]
    
    pointM0 =  int(ground.split()[0])
    pointN0 =  int(ground.split()[1])
    calc0 =  int(ground.split()[2])
    E0 = ground.split()[3]
    G0 = ground.split()[4]
    reZ10 = ground.split()[5]
    imZ10 = ground.split()[6]

    
    pointM =  int(upper.split()[0])
    pointN =  int(upper.split()[1])
    calc =  int(upper.split()[2])
    E = upper.split()[3]
    G = upper.split()[4]
    reZ1 = upper.split()[5]
    imZ1 = upper.split()[6]
    
    
    if(mpf(F) == 0):    
        z0 = mpf("1000")
    else:    
        z0 = 1/sqrt(mpf(F))
            
    maxM = min(pointM ,pointM0)    
    
    
    
    
    if((abs(m)==1) or (m==0)):

        p0 = WaveFunc(0,0,0,pointM0,pointN0,E0,G0,reZ10,imZ10,F,num0)
        num0 = p0.mode

        p = WaveFunc(n1,n2,m,pointM,pointN,E,G,reZ1,imZ1,F,num)
        num = p.mode
        
        MN0 = p0.MN
        MN = p.MN

        if(abs(m)==1):
            tr_x = pi*quad(lambda mu: quad(lambda nu:   conj(MN0(mu,nu))*(mu*nu)*MN(mu,nu)     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM])
            tr_y = m*mpc(j)*tr_x
            tr_z = 0

        if(m==0):
            tr_x = 0
            tr_y = 0
            tr_z = 2*pi*quad(lambda mu: quad(lambda nu:   conj(MN0(mu,nu))*((mu*mu - nu*nu)/2)*MN(mu,nu)     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM])            

    else:
        tr_x = 0
        tr_y = 0
        tr_z = 0

    
    f = open('000--->%i%i%i.txt'%(n1,n2,m),'a')
    print >>f, F.ljust(30),   nstr(tr_x,20).ljust(60),  nstr(tr_y,20).ljust(60),  nstr(tr_z,20).ljust(60),   str(num0)+str(num)
    f.close()
        
    
        
        
