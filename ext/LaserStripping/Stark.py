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
        self.array_E = []
        self.array_G = []
        self.pass_calc = True
        self.step_E = 0
        self.step_G = 0
        self.tempE = mpf("1e100")
        self.tempG = mpf("1e100")




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
        
        self.tempE = mpf("1e100")
        self.tempG = mpf("1e100")
        
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
        
        if (F == mpf("0.0")):
            self.Gamma = mpf("0.0")
        else:
            R = power(-2*self.Energy,fdiv(3,2))/F
            self.Gamma = power(4*R,2*n2+m+1)/(n*n*n*fac(n2)*fac(n2+m))*exp(-2*R/3-(n*n*n*F/4)*(34*n2*n2+34*n2*m+46*n2+7*m*m+23*m+fdiv(53,3)))
        
        self.step_E = fabs((self.Energy + fdiv(1,2)/(n*n))/100)
        self.step_G = fabs(self.Gamma/100)
       
        return
    
    
    
    
    
    def predict_EG(self):
        
        k = 100
        
        n = len(self.array_E)
        m = min(n,k)
        self.Energy = 0
        for i in range(0,m):
            self.Energy += power(-1,m+i+1)*fac(m)/(fac(m-i)*fac(i))*self.array_E[i+n-m]


        n = len(self.array_G)
        m = min(n,k)
        sumG = 0
        for i in range(0,m):
            sumG += power(-1,m+i+1)*fac(m)/(fac(m-i)*fac(i))*log(self.array_G[i+n-m])
            
        self.Gamma = exp(sumG)

        return
             
             
             
    
    def initialEG(self):
    
        
        if(len(self.array_G)==1):
            if(self.array_G[0]==mpf("0.0")):
                self.array_G.pop(0)
        
        if(len(self.array_G)>1):
            if(self.array_G[1]/self.array_G[0] > mpf("1e10")):
                self.array_G.pop(0)
                
        
        if(self.pass_calc or (len(self.array_E)<5) or (len(self.array_G)<5)):
            self.def_start_EG()
            print "Epred_an = ", self.Energy, "  Gpred_an = ", self.Gamma, "  dE_an = ",self.step_E,"  dG_an = ",self.step_G

        else:

            real_last_E = self.array_E[len(self.array_E)-1]
            real_last_G = self.array_G[len(self.array_G)-1]
            
            self.array_E.pop(len(self.array_E)-1)
            self.array_G.pop(len(self.array_G)-1)
            
            self.predict_EG()
            
            predicted_last_E = self.Energy
            predicted_last_G = self.Gamma
            
            self.array_E.append(real_last_E)
            self.array_G.append(real_last_G)
            
            self.predict_EG()
            
            next_predicted_E = self.Energy
            next_predicted_G = self.Gamma         
                        
            self.step_E = fabs(next_predicted_E*(real_last_E - predicted_last_E)/real_last_E)/10
            self.step_G = fabs(next_predicted_G*(real_last_G - predicted_last_G)/real_last_G)/10    
            
            print "Epred = ", self.Energy, "  Gpred = ", self.Gamma, "  dE = ",self.step_E,"  dG = ",self.step_G
    
        return





    def find_Z1(self):
        
    
        s = Simplex(self.absM, [re(self.Z1), im(self.Z1)], [mpf("1e-5"),mpf("1e-5")])
        (values, err, iter) = s.minimize(power(10,self.err_exp), 10000000,0)
        self.Z1 = mpc(values[0], values[1])
        return 

    
      
        
    def absB(self,args):
        
        err_GE = mpf("1e-20")
        
        self.Energy = args[0]
        self.Gamma = args[1]
        
        
        self.find_Z1()
        self.Z2 = fsub(4,self.Z1)
        par = self.a.getB(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
        line =  par.replace("(","").replace(")","").rsplit(" ")
        B = mpc(line[0], line[1])
#        print  "E = ",nstr(self.Energy,30),"   G = ",nstr(self.Gamma,30),"  B = ",nstr(B,30),"  absB = ",nstr(fabs(B),30)
        
        if((fabs((self.tempE - self.Energy)/self.Energy) < err_GE) and (fabs((self.tempG - self.Gamma)/self.Gamma) < err_GE)):
            return mpf("1e-1000000")
        else:
#            print fabs((self.tempE - self.Energy)/self.Energy)
#            print fabs((self.tempG - self.Gamma)/self.Gamma)
            
            self.tempE = self.Energy
            self.tempG = self.Gamma

            
            return fabs(B)
    
    
    
    
    def find_EG(self):
        
        Ein = self.Energy
        Gin = self.Gamma
        
        s = Simplex(self.absB, [self.Energy, self.Gamma], [self.step_E,self.step_G])
        (values, err, iter) = s.minimize(mpf("1e-10000"), 10000000,0)

        self.pass_calc = (err==mpf("0.0"))
        
        if (self.pass_calc):
            self.Energy = Ein
            self.Gamma = Gin
            
        self.array_E.append(self.Energy)
        self.array_G.append(self.Gamma)
        
        return 

    
        
            
#E =  -0.0609859323708186145609132465066    G =  0.00000201579857103347930193967475202   B =  (5.06916915324052272960434397358e-37 - 3.31449855012918010253005911196e-37j)   absB =  6.05659776962063704430148146337e-37
    



n1 = 0
n2 = 0
m = 0
point_gs = 10
err_exp = -50








b = Stark_calc(n1, n2, m, point_gs,err_exp)



for i in range(1,300):
    b.F = i*mpf("1.0e-5")
    b.defr_parameters_forF()
    b.initialEG()
    b.find_EG()

#    b.F = 10*mpf("1.0e-5")
#    mp.prec = b.a.calcPrecisionForM(str(b.F),"("+str(b.Energy)+" "+ str(b.Gamma*mpf("-0.5"))+")", "("+str(re(b.Z1))+" "+ str(im(b.Z1))+")")
#    b.find_Z1()
#    self.Z2 = fsub(4,b.Z1)
#    self.a.calcPrecisionForN(str(self.F),"("+str(self.Energy)+" "+ str(self.Gamma*mpf("-0.5"))+")", "("+str(re(self.Z2))+" "+ str(im(self.Z2))+")")
        
    
    print "F = ", b.F, "  Energy = ",b.Energy,"  G = ",b.Gamma,"  pass_calc = ", b.pass_calc


"""
d = mpf("0.001")

k = 100
for i in range(0,k):
    b.F = d*i
    b.initialEG()
    b.array_E.append(b.Energy)
    b.array_G.append(b.Gamma)
    


b.F = d*k
b.initialEG()

#print "real = ",b.Energy
print "real = ",b.Gamma

b.pass_calc = False
b.initialEG()
    
#print "pred = ",b.Energy
print "pred = ",b.Gamma
print "err  = ", b.step_G

#for i in range(0,k-5):
#    print b.array_G[i]
"""



print "mp.dps =",mp.dps


#print b.def_point1(10)



#(2.000001999995500022249843345128907569821710508245174093240455103836507433670436765193680160603367165824627094388

