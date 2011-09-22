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







(n1,n2,m) = (0,3,1)

num0 = 0
num = 0


for k in xrange(0,1000000000):
    if not os.path.exists("_run"):    sys.exit(0)

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
        num0 += p0.switch

        p = WaveFunc(n1,n2,m,pointM,pointN,E,G,reZ1,imZ1,F,num)
        num += p.switch
        
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
        
    if(m==0):
        orth = 2*pi*quad(lambda mu: quad(lambda nu:   conj(MN0(mu,nu))*MN(mu,nu)     *nu*mu*(nu*nu+mu*mu),[0,sqrt(mu*mu+2*z0)]),[0,maxM])
    else:
        orth = 0

    
    f = open('000--->%i%i%i.txt'%(n1,n2,m),'a')
    print >>f, F.ljust(30),   nstr(tr_x,20).ljust(60),  nstr(tr_y,20).ljust(60),  nstr(tr_z,20).ljust(60)
    f.close()
        
    ch = open('ch%i%i%i.txt'%(n1,n2,m),'a')
    print >>ch, F.ljust(30),  (str(num0)+str(num)).ljust(10),   nstr(orth,20).ljust(60)  
    ch.close()
    

        
        