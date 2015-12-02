import sys
#import math
import random
import numpy
import argparse
import walker
import copy
import csv


# Some parsers for command line options
parser = argparse.ArgumentParser(description='Walk with killing.')
parser.add_argument('infofile',  metavar='infofile',  nargs='?',  help='output file 1',  default='info')
parser.add_argument('distfile',  metavar='distfile',  nargs='?',  help='output file 2',  default='dist')
parser.add_argument('-T', metavar='T',  nargs='?',  type=int,  help='total runtime', default=10)
parser.add_argument('-s', metavar='s',  nargs='?', type=int,  help='number of discretization steps', default=100)
parser.add_argument('-w', metavar='w',  nargs='?', type=int,  help='total number of walkers', default=10)
parser.add_argument('-d', metavar='d',  nargs='?', type=int,  help='dimension of hypercube graph', default=10)
args=parser.parse_args()

#script, filename1, filename2 = argv
#print("usage: python adiabaticv2.py info_file distribution_file")
infofile = open(args.infofile, 'w') # writing file
infofile.truncate()
writer = csv.writer(infofile)
distfile = open(args.distfile, 'w')  # writing file
distfile.truncate()
##################### Parameters
n=args.d # dimension of hypercube
T=(n**2)
deltaT=1/float(3)
timesteps=int(T/deltaT) # number of steps
if (T/timesteps > 1): 
    print("badness")
    sys.exit()

vertexnumber=int(2**n) # for hypercube
initialwalkernum=n**2 # initial number of walkers should be logarithmic in the number of vertices
##################### Use NetworkX to generate a hypercube; generate various dictionaries with the output
       
def revive( walkerList ):
    tmpList = list()
    tmpList.extend(walkerList)
    while len(tmpList) < initialwalkernum:
        w = random.choice(walkerList)                           # choose a random surviving member
        tmpList.append(w.spawn())                            # revive one
    return tmpList
    
binList = list(range(n+1))
def histogram( walkerList ):
    histoList = list()
    for w in walkerList:
        histoList.append(w.potential())
    return numpy.histogram(histoList, bins=binList, density=False)
    
survivors=list()
for j in range(initialwalkernum):                                       # randomize initial distribution
    rando = random.randint(0,vertexnumber-1)
    survivors.append(walker.Walker(rando, n))                 #add a new walker
    
percentDone = 0.000
for t in range(timesteps): # random walk loop
    s=t/float(timesteps) # adiabatic scheduling parameter
    
    writer.writerow((histogram(survivors)[0]))
    distfile.write(" %s \n" % (histogram(survivors)[0]))
    
    if s > percentDone:
        sys.stdout.write("\r %s %% Complete" % (percentDone*100))
        percentDone+=0.001
        sys.stdout.flush()
    
    walkers=survivors
    survivors=list() #Note: this comes after the revival since old/new keep revivals separate
    if len(walkers)<initialwalkernum:
        walkers = revive(walkers)
    
    
    for v in walkers: # move walkers, one walker at a time
            
        probedge= (1-s)*deltaT
        probdead= v.potential()*deltaT*s
        probstay = 1- probdead - probedge
        
        randprob=random.random()        
        
        if(randprob<=probedge): # walk to a neighbor
            v.walk()
            survivors.append(v)
        elif(randprob <= probedge+probstay): # do nothing
            survivors.append(v)
        # if you weren't caught in the last two ifs, you die or, in other words, we just don't add you as a survivor
        
    if len(survivors) == 0:
        print("All walkers died")
        break
        

flag = False
for v in survivors:
    if v.vertex==0:
        flag=True
        break

if flag:
    print("\n Solution Found")
else:
    print("\n Failed to find solution")
