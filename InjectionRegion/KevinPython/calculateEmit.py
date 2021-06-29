
from orbit.teapot.teapot import NodeTEAPOT
import math
from bunch import BunchTwissAnalysis
class Calc_Emit(NodeTEAPOT):
    """
    Generates the Histograms for distribution in x,y, and x=y directions.
    The z direction from the original WS_Node is repurposed for the diagonal x=y.
    It is not parallel !!!!
    """
    def __init__(self,name = "CalcEmit",toFile=False,fileName="output.txt"):
        NodeTEAPOT.__init__(self,name)
        self.toFile=toFile
        self.fileName=fileName

    def track(self, paramsDict):
        """
        The WS_AccNode class implementation of the AccNode class track(probe) method.
        """
        debug=False
        bunch = paramsDict["bunch"]
        gamma = bunch.getSyncParticle().gamma()
        beta = bunch.getSyncParticle().beta()
        
        tot_x=0
        tot_px=0
        tot_y=0
        tot_py=0
        
        tot_x2=0
        tot_px2=0
        tot_x_px=0
        
        tot_y2=0
        tot_py2=0
        tot_y_py=0       
        if bunch.getSize() >1:
		for i in range(bunch.getSize()):   
			tot_x=tot_x+bunch.x(i)*1000.
			tot_x2=tot_x2+bunch.x(i)*1000.*bunch.x(i)*1000.
			tot_px=tot_px+bunch.px(i)*1000.
			tot_px2=tot_px2+bunch.px(i)*1000.*bunch.px(i)*1000.
			tot_y=tot_y+bunch.y(i)*1000.
			tot_y2=tot_y2+bunch.y(i)*1000.*bunch.y(i)*1000.
			tot_py=tot_py+bunch.py(i)*1000.
			tot_py2=tot_py2+bunch.py(i)*1000.*bunch.py(i)*1000.
			tot_x_px=tot_x_px+bunch.x(i)*1000.*bunch.px(i)*1000.
			tot_y_py=tot_y_py+bunch.y(i)*1000.*bunch.py(i)*1000.       	
		avg_x=tot_x/bunch.getSize()
		avg_px=tot_px/bunch.getSize()
		avg_y=tot_y/bunch.getSize()
		avg_py=tot_py/bunch.getSize() 
		
		avg_x2=tot_x2/bunch.getSize()
		avg_px2=tot_px2/bunch.getSize()
		avg_x_px=tot_x_px/bunch.getSize()
		
		avg_y2=tot_y2/bunch.getSize()
		avg_py2=tot_py2/bunch.getSize()
		avg_y_py=tot_y_py/bunch.getSize()  
		
		var_x=avg_x2-math.pow(avg_x,2)
		var_px=avg_px2-math.pow(avg_px,2)
		var_x_px=avg_x_px-avg_x*avg_px
		
		var_y=avg_y2-math.pow(avg_y,2)
		var_py=avg_py2-math.pow(avg_py,2)
		var_y_py=avg_y_py-avg_y*avg_py        
		if debug==True:
			print "(tot_x2) = %f"%(tot_x2)
			print "(tot_px2) = %f"%(tot_px2)
			print "(tot_x_px) = %f"%(tot_x_px)
			print "(avg_x2) = %f"%(avg_x2)
			print "(avg_px2) = %f"%(avg_px2)
			print "(avg_x_px) = %f"%(avg_x_px)        	
			print "(emit x)^2 = %f"%(avg_x2*avg_px2-avg_x_px*avg_x_px)
			print "(tot_y2) = %f"%(tot_y2)
			print "(tot_py2) = %f"%(tot_py2)
			print "(tot_y_py) = %f"%(tot_y_py)
			print "(avg_y2) = %f"%(avg_y2)
			print "(avg_py2) = %f"%(avg_py2)
			print "(avg_y_py) = %f"%(avg_y_py)        	
			print "(emit y)^2 = %f"%(avg_y2*avg_py2-avg_y_py*avg_y_py)    
			print "(var_x) = %f"%(var_x)
			print "(var_px) = %f"%(var_px)
			print "(var_x_px) = %f"%(var_x_px)
			print "(var_y) = %f"%(var_y)
			print "(var_py) = %f"%(var_py)
			print "(var_y_py) = %f"%(var_y_py)
			print "(var_x*var_px-var_x_px*var_x_px) = %f"%(var_x*var_px-var_x_px*var_x_px)   
		emit_x=math.sqrt(var_x*var_px-var_x_px*var_x_px)
		emit_y=math.sqrt(var_y*var_py-var_y_py*var_y_py)
		
		twiss_analysis = BunchTwissAnalysis()        
		node = paramsDict["node"]
		bunch = paramsDict["bunch"]
		pos = paramsDict["path_length"]
	
		gamma = bunch.getSyncParticle().gamma()
		beta = bunch.getSyncParticle().beta()
		twiss_analysis.analyzeBunch(bunch)
		x_rms = math.sqrt(twiss_analysis.getTwiss(0)[1]*twiss_analysis.getTwiss(0)[3])*1000.
		y_rms = math.sqrt(twiss_analysis.getTwiss(1)[1]*twiss_analysis.getTwiss(1)[3])*1000.
		z_rms = math.sqrt(twiss_analysis.getTwiss(2)[1]*twiss_analysis.getTwiss(2)[3])*1000.
		xp_rms= math.sqrt(twiss_analysis.getTwiss(0)[2]*twiss_analysis.getTwiss(0)[3])*1000.
		yp_rms= math.sqrt(twiss_analysis.getTwiss(1)[2]*twiss_analysis.getTwiss(1)[3])*1000.
		nParts = bunch.getSizeGlobal()
		(alphaX,betaX,emittX) = (twiss_analysis.getTwiss(0)[0],twiss_analysis.getTwiss(0)[1],twiss_analysis.getTwiss(0)[3]*1.0e+6)
		(alphaY,betaY,emittY) = (twiss_analysis.getTwiss(1)[0],twiss_analysis.getTwiss(1)[1],twiss_analysis.getTwiss(1)[3]*1.0e+6)
		norm_emittX = emittX*gamma*beta
		norm_emittY = emittY*gamma*beta
		#---- phi_de_emittZ will be in [pi*deg*MeV]
		eKin = bunch.getSyncParticle().kinEnergy()*1.0e+3
	
		s_prt = " %35s  %4.5f "%(node.getName(),pos)
		s_prt += "  %5.3f  %5.3f  "%(x_rms,y_rms)
		s_prt += "  %10.6f   %8d "%(eKin,nParts)
		print s_prt	        
		if (self.toFile==False):
		 
			print " (emit x, emit y)= (%f,%f) " %(emit_x,emit_y)
		else:
		    fileOut=open(self.fileName,'a')
		    #print self.getName()
		    s = " %35s  %4.5f \n"%(node.getName(),pos)
		    s += "   %6.4f  %6.4f  %6.4f  %6.4f   \n"%(alphaX,betaX,emittX,norm_emittX)
		    s += "   %6.4f  %6.4f  %6.4f  %6.4f   \n"%(alphaY,betaY,emittY,norm_emittY)
		    s += "   %5.3f  %5.3f \n"%(x_rms,y_rms)
		    s += "   %5.3f  %5.3f \n"%(xp_rms,yp_rms)
		    s += "  %10.6f   %8d "%(eKin,nParts)
		    fileOut.write(s +"\n")
		    fileOut.flush()            
		    fileOut.write(" (x rms, y rms)= (%f,%f) \n" %(math.sqrt(var_x),math.sqrt(var_y)))
		    fileOut.write(" (xp rms, yp rms)= (%f,%f) \n" %(math.sqrt(var_px),math.sqrt(var_py)))
		    fileOut.write(" (x avg, y avy)= (%f,%f) \n" %(avg_x,avg_y))
		    fileOut.write(" (xp avg, yp avg)= (%f,%f) \n" %(avg_px,avg_py))		    
		    fileOut.write(" (emit x, emit y)= (%f,%f) \n" %(emit_x,emit_y))
		    fileOut.write(" (normal emit x, normal emit y)= (%f,%f) \n" %(emit_x*beta*gamma,emit_y*beta*gamma))
		    fileOut.close()
	elif bunch.getSize() ==1:
	    fileOut=open(self.fileName,'a')
	    #print self.getName()
	    node = paramsDict["node"]
	    pos = paramsDict["path_length"]
	    s = " %35s  %4.5f \n"%(node.getName(),pos)
	    fileOut.write(s +"\n")
	    fileOut.flush()       
	    fileOut.write("pencilbeam \n")
	    fileOut.write(" (x rms, y rms)= (%f,%f) \n" %(0,0))
	    fileOut.write(" (xp rms, yp rms)= (%f,%f) \n" %(0,0))
	    fileOut.write(" (x avg, y avg)= (%f,%f) \n" %(bunch.x(0),bunch.y(0)))
	    fileOut.write(" (xp avg, yp avg)= (%f,%f) \n" %(bunch.px(0),bunch.py(0)))		    
	    fileOut.write(" (emit x, emit y)= (%f,%f) \n" %(0,0))
	    fileOut.write(" (normal emit x, normal emit y)= (%f,%f) \n" %(0,0))
	    fileOut.close()		