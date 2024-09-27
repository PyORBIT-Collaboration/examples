from ROOT import *
import sys
import os, argparse
import signal

def signal_handler(signal, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)

parser = argparse.ArgumentParser(description="%prog [options]", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument("--directory", dest='directory', default="/share/t3data3/khildebr/RandomConeV204MC/", help="Path to the input directory")
parser.add_argument("--fileName", dest='fileName', default="output5.txt", help="Input File to read")
#parser.add_argument("--totalTGraph", type=int, dest='totalTGraph', default=-1, help="the number of tgraphs to include in total tgraph")
parser.add_argument("--totalTGraph", type=int, dest='totalTGraph', default=0, help="the number of tgraphs to include in total tgraph")
parser.add_argument("--xmin", type=float, dest='xmin', default=-.065, help="xmin for tgraph")
parser.add_argument("--xmax", type=float, dest='xmax', default=.065, help="xmax for tgraph")
parser.add_argument("--ymin", type=float, dest='ymin', default=-.065, help="ymin for tgraph")
parser.add_argument("--ymax", type=float, dest='ymax', default=.065, help="ymax for tgraph")
parser.add_argument("--xaxis", type=int, dest='xaxis', default=0, help="what to plot on xaxis for tgraph. Options are [0-6]=(x,px,y,py,z,pz,s)")
parser.add_argument("--yaxis", type=int, dest='yaxis', default=2, help="what to plot on yaxis for tgraph. Options are [0-6]=(x,px,y,py,z,pz,s)")
parser.add_argument("--directory", dest='directory', default="output5", help="directory to put graphs into")
parser.add_argument("--imageType", dest='imageType', default="png", help="what file type to save image as (ie png)")
#parser.add_argument("--markerSize",type=float, dest='markerSize', default="1", help="size of markers in graph")
parser.add_argument("--markerSize",type=float, dest='markerSize', default=".25", help="size of markers in graph")
parser.add_argument("--markerStyle",type=int, dest='markerStyle', default="8", help="style of markers in graph")
args = parser.parse_args()

doAll=True
numberToDo=-1
if (args.totalTGraph>0):
    doAll=False
    numberToDo=args.totalTGraph
    
#fileName=sys.argv[1]
if not os.path.isdir(args.directory):
    os.mkdir(args.directory)
    #print("%s directory does not exist"%(args.directory))
    #sys.exit(0)
theNaming=["X","PX","Y","PY","Z","PZ","S"]
theXAxisVariable=theNaming[args.xaxis]
theYAxisVariable=theNaming[args.yaxis]
if not os.path.isdir("%s/%svs%s"%(args.directory,theXAxisVariable,theYAxisVariable)):
    os.mkdir(("%s/%svs%s"%(args.directory,theXAxisVariable,theYAxisVariable)))
openedFile=open(args.fileName,'r')

#stops each canvas from being physically drawn
gROOT.SetBatch(True)

lines=openedFile.readlines()
graphs=[]
oneGraph= TGraph()
totalGraph=TGraph()
#count number of tgraphs
counter=0
#count number of lines read
counterLinesRead=0
for line in lines:
    #print line[0]
    #print line
    tokens=line.split("(")
    #print tokens[0]
    #print tokens[1]
    #print tokens[2]
    if (int(tokens[0].strip())==0 and counterLinesRead!=0):
        theCanvas=TCanvas("TGraph","TGraph",0,0,500,500)
        oneGraph.GetYaxis().SetRangeUser(args.ymin,args.ymax);
        oneGraph.GetXaxis().SetLimits(args.xmin,args.xmax);
        oneGraph.SetMarkerStyle(args.markerStyle)
        oneGraph.SetMarkerSize(args.markerSize)
        oneGraph.SetTitle("theTitle")
        oneGraph.Draw("AP")   
        theCanvas.Print("%s/%svs%s/temp%d.%s"%(args.directory,theXAxisVariable,theYAxisVariable,counter,args.imageType))
        #theCanvas.SaveAs("temp%d.png"%counter)
        theCanvas.Clear()
        oneGraph= TGraph()
        counter+=1
    counterLinesRead+=1
    #tokens[2]=tokens[2].strip(')')
    tokens[2]=tokens[2].strip().strip(')')
    #parameters[0-5]=(x,px,y,py,z,pz)
    #parameters[0]=x
    #parameters[1]=px
    #parameters[2]=y
    #parameters[3]=py
    #parameters[4]=z
    #parameters[5]=pz
    parameters=tokens[2].split(',')
    #print tokens[2]
    oneGraph.SetPoint(oneGraph.GetN(),float(parameters[args.xaxis]),float(parameters[args.yaxis]))
    if (doAll==True or numberToDo>counter):
        totalGraph.SetPoint(totalGraph.GetN(),float(parameters[args.xaxis]),float(parameters[args.yaxis]))

#still need to print last graph
theCanvas=TCanvas("TGraph","TGraph",0,0,500,500)
oneGraph.GetYaxis().SetRangeUser(args.ymin,args.ymax);
oneGraph.GetXaxis().SetLimits(args.xmin,args.xmax);
oneGraph.SetMarkerStyle(args.markerStyle)
oneGraph.SetMarkerSize(args.markerSize)
oneGraph.SetTitle("theTitle")
oneGraph.Draw("AP")   
theCanvas.Print("%s/%svs%s/temp%d.%s"%(args.directory,theXAxisVariable,theYAxisVariable,counter,args.imageType))
#theCanvas.SaveAs("temp%d.png"%counter)
theCanvas.Clear()
oneGraph= TGraph()
counter+=1

#save totalgraph
theCanvas=TCanvas("TGraph","TGraph",0,0,500,500)
totalGraph.GetYaxis().SetRangeUser(args.ymin,args.ymax);
totalGraph.GetXaxis().SetLimits(args.xmin,args.xmax);
totalGraph.SetMarkerStyle(args.markerStyle)
totalGraph.SetMarkerSize(args.markerSize)
totalGraph.SetTitle("theTitle")
totalGraph.Draw("AP")   
theCanvas.Print("%s/%svs%s/aTotal%d.%s"%(args.directory,theXAxisVariable,theYAxisVariable,counter,args.imageType))
#theCanvas.SaveAs("temp%d.png"%counter)
theCanvas.Clear()
totalGraph= TGraph()
