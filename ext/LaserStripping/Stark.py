from starkeffect import *
from ext.las_str.SimplexModMPF import Simplex
from mpmath import *









class Stark_calc: 
    
    

    def __init__(self,_n1, _n2, _m, point_ground_state, _err_exp):
        
        mp.prec = 10000
        self.n1 = _n1
        self.n2 = _n2
        self.m = abs(_m)
        self.n = _n1 + _n2 + abs(_m) + 1
        self.Energy = fdiv(-1,2*self.n*self.n)
        self.Gamma = mpf("0")
        self.Z1 = mpc(fdiv(2*(2*n1 + abs(m) + 1), self.n), 0)
        self.Z2 = fsub(4,self.Z1)
        self.F = mpf(0)
        self.a = Functions(n1, n2, m, self.def_point(point_ground_state),_err_exp)
        self.err_exp = _err_exp




    def def_point(self, p1_for_ground):

        
        level0 = exp(p1_for_ground*p1_for_ground*mpf("-0.5"))*sqrt(p1_for_ground)
        for i in range(100,0,-1):            
            level = exp(i*i*mpf("-0.5")/self.n)*sqrt(i)*power(i,self.m)*hyp1f1(-self.n1, 1 + self.m, fdiv(i*i,self.n))
            if (abs(level)>level0):
                break
        return i + 1


    def absM(self,args):
        
        self.Z1 = mpc(args[0],args[1])
        par = self.a.getM(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z1))+" "+ str(im(self.Z1))+")")
        line =  par.replace("(","").replace(")","").rsplit(" ")
        comp = mpc(line[0], line[1])
#        print "args= ",nstr(re(self.Z1),20),nstr(im(self.Z1),20),"copm= ",nstr(fabs(comp),20)
        return fabs(comp)
    
    
    def defr_parameters_forF(self):
        
        mp.prec = self.a.calcPrecisionForM(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z1))+" "+ str(im(self.Z1))+")")
        self.find_Z1()
        self.Z2 = fsub(4,b.Z1)
        self.a.calcPrecisionForN(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
        
        return


    
    def def_start_EG(self):
        
        n1 = self.n1
        n2 = self.n2
        m = self.m
        n = self.n
        F = self.F
        
        self.Energy = (
        -fdiv(1,2)/(n*n)
        +F*fdiv(3,2)*n*(n1-n2)
        -F*F*fdiv(1,16)*power(n,4)*(17*n*n-3*(n1-n2)*(n1-n2)-9*m*m+19)
        +F*F*F*fdiv(3,32)*power(n,7)*(n1-n2)*(23*n*n-(n1-n2)*(n1-n2)+11*m*m+39)
        -F*F*F*F*fdiv(1,1024)*power(n,10)*(5487*power(n,4)+35182*n*n+(5754-1134*m*m+1806*n*n)*(n1-n2)*(n1-n2)-549*m*m*m*m-3402*n*n*m*m+147*power(n1-n2,4)-8622*m*m+16211)
        +F*F*F*F*F*fdiv(3,1024)*power(n,13)*(n1-n2)*(10563*n*n*n*n+90708*n*n+(n1-n2)*(n1-n2)*(220*m*m+98*n*n+780)+772*n*n*m*m-21*power(n1-n2,4)+725*m*m*m*m+830*m*m+59293)
        )
        
        R = power(-2*self.Energy,fdiv(3,2))/F
        
        self.Gamma = power(4*R,2*n2+m+1)/(n*n*n*fac(n2)*fac(n2+m))*exp(-2*R/3-(n*n*n*F/4)*(34*n2*n2+34*n2*m+46*n2+7*m*m+23*m+fdiv(53,3)))
       
        return
             
    

    def find_Z1(self):
        
    
        s = Simplex(self.absM, [re(self.Z1), im(self.Z1)], [mpf("1e-5"),mpf("1e-5")])
        (values, err, iter) = s.minimize(power(10,self.err_exp), 10000000,0)
        self.Z1 = mpc(values[0], values[1])
        return 

    
      
        
    def absB(self,args):
        
        self.Energy = args[0]
        self.Gamma = args[1]
        
        
        self.find_Z1()
        self.Z2 = fsub(4,self.Z1)
        par = self.a.getB(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
        line =  par.replace("(","").replace(")","").rsplit(" ")
        B = mpc(line[0], line[1])
        print  "E = ",nstr(self.Energy,30),"   G = ",nstr(self.Gamma,30),"  B = ",nstr(B,30),"  absB = ",nstr(fabs(B),30)
        
        return fabs(B)
    
    
    
    
    def find_EG(self):
        
        s = Simplex(self.absB, [self.Energy, self.Gamma], [mpf("1e-2"),mpf("1e-2")])
        (values, err, iter) = s.minimize(mpf("1e-50"), 10000000,0)

        return 

    
        
            

    



n1 = 0
n2 = 1
m = 1
point_gs = 10
err_exp = 0










b = Stark_calc(n1, n2, m, point_gs,err_exp)


b.F = mpf("5.5e-3")
b.def_start_EG()
b.defr_parameters_forF()
b.find_EG()






print "mp.dps =",mp.dps


#print b.def_point1(10)



#(2.000001999995500022249843345128907569821710508245174093240455103836507433670436765193680160603367165824627094388

