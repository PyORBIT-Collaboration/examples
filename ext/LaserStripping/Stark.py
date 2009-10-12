from starkeffect import *
from ext.las_str.SimplexModMPF import Simplex
from mpmath import *









class Stark_calc: 
    
    

    def __init__(self,_n1, _n2, _m, _point1, _err_exp):
        
        mp.prec = 10000
        self.n1 = _n1
        self.n2 = _n2
        self.m = abs(_m)
        self.n = _n1 + _n2 + abs(_m) + 1
        self.Energy = fdiv(-1,2*self.n*self.n)
        self.Gamma = mpf(0)
        self.Z1 = mpc(fdiv(2*(2*n1 + abs(m) + 1), self.n), 0)
        self.Z2 = fsub(4,self.Z1)
        self.F = mpf(0)
        self.a = Functions(n1, n2, m, point1,_err_exp)
        self.err_exp = _err_exp




    def def_point1(self, p1_for_ground):

        
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
       

    

    def find_Z1(self):
        
    
        s = Simplex(self.absM, [re(self.Z1), im(self.Z1)], [mpf("1e-5"),mpf("1e-5")])
        (values, err, iter) = s.minimize(power(10,self.err_exp), 10000000,0)
        self.Z1 = mpc(values[0], values[1])
        return 

    
      
        
    def absB(self,args):
        
        self.Energy = args[0]
        self.Gamma = args[1]
        
#        mp.prec = self.a.calcPrecisionForM(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z1))+" "+ str(im(self.Z1))+")")
#        self.find_Z1()
        self.Z1 = mpc(fdiv(2*(2*n1 + abs(m) + 1), self.n), 0)
        self.Z2 = fsub(4,self.Z1)
        print "flag"
        N_prec, N, derN = self.a.getN(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
#        print "mp.prec, N_prec = ",mp.prec, N_prec
#        a_prec, a, der_a = self.a.get_a(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
#        print "mp.prec, a_prec = ",mp.prec, a_prec
#        b_prec, b, der_b = self.a.get_b(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
#        print "mp.prec, b_prec = ",mp.prec, b_prec
        
        print N, derN
#        print a, der_a
#        print b, der_b
        
        return 
        
            

    



n1 = 0
n2 = 0
m = 0
point1 = 20
err_exp = 0







mp.prec = 10000




b = Stark_calc(n1, n2, m, point1,err_exp)
b.F = mpf("1e-4")




b.absB([b.Energy, b.Gamma])


print "mp.dps =",mp.dps


#print b.def_point1(10)



#(2.000001999995500022249843345128907569821710508245174093240455103836507433670436765193680160603367165824627094388

