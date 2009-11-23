from ext.las_str.StarkEffect  import *
from mpmath import *




(n1,n2,m) = (0,0,0)


b = Stark_calc(n1, n2, m, 10)



for i in xrange(100,10000000):
    b.F = i*mpf("1.0e-6")
    b.main()


    
    f = open('%i%i%i.txt'%(n1,n2,m),'a')
    print >>f, str(b.F),"\t",nstr(b.Energy,20),"\t",nstr(b.Gamma,20)
    f.close()
    
    f = open('_%i%i%i.txt'%(n1,n2,m),'a')
    print >>f, b.pointM,"\t",b.pointN,"\t",b.calc,"\t",str(b.Energy),"\t",str(b.Gamma),"\t",str(re(b.Z1)),"\t",str(im(b.Z1)),"\t",nstr(b.F,50)
    f.close()








