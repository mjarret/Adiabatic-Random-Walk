import sys
#import math
import random
import numpy
import networkx as nx
import time
import argparse
import walker
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
n=17 # dimension of hypercube
T=1000*(n**2)
deltaT=1
timesteps=int(T/deltaT) # number of steps
if (T/timesteps > 1): 
    print("badness")
    sys.exit()

vertexnumber=2**n # for hypercube
initialwalkernum=n**2 # initial number of walkers should be logarithmic in the number of vertices
##################### Use NetworkX to generate a hypercube; generate various dictionaries with the output
print("generating graph")
hyper=nx.hypercube_graph(n) # generate hypercube using NetworkX
print("done")

vertexlist=hyper.nodes()
edgelist=hyper.edges()
adjlist=hyper.adjacency_list()
incidencedict=dict() # adjacency dictionary
hammingdist=dict() # hamming weight dictionary
layerdict=dict() # dictionary of vertices organized according to hamming weights
#eg for n=3 layerdict[2]={(1,1,0),(1,0,1),(0,1,1)}

# Function will be used to repeatedly calculate Hamming weight
def hammingWeight( vertex ):
    c=0
    for v in vertex:
        c += v
    return c
    
def potential( vertex ):    
    return hammingWeight(vertex)/float(n)
    
oneUp = list()
for i in range(n):
    oneUp.append(2**i)
print(oneUp)

def walk( vertex ):
    vertex += random.choice()
    
def revive( walkerList ):
    tmpList = list()
    tmpList.extend(walkerList)
    while len(tmpList) < initialwalkernum:
        w = random.choice(walkerList)    # choose a random surviving member
        tmpList.append(copy.copy(w))                       # revive one
    return tmpList
    
binList = list(range(n+1))
def histogram( walkerList ):
    histoList = list()
    for v in walkerList:
        histoList.append(hammingWeight(v))
    return numpy.histogram(histoList, bins=binList, density=True)

for dist in range((2**n)+1):
    layerdict[dist]=list()

for v in range(vertexnumber):
    vert=vertexlist[v]
    adj=adjlist[v]
    incidencedict[vert]=adj
    weight=0 # variable for hamming weight of v
    for ind in range(n):
        weight=weight+vert[ind]
    hammingdist[vert]=weight
    layerdict[weight].append(vert)
##################### MCMC Walk
#walkerlistold=list() # list of walker locations
survivors=list()

for j in range(initialwalkernum): # randomize initial distribution
    rando = random.randint(0,vertexnumber-1)
    survivors.append(vertexlist[rando])
thewalkingdead=0 # initialize number of walkers that died in the last step
time1=time.time() # for timing of MCMC

percentDone = 0.000
for t in range(timesteps): # random walk loop
    s=t/float(timesteps) # adiabatic scheduling parameter
    if s > percentDone:
        sys.stdout.write("\r %s %% Complete" % percentDone*100)
        percentDone+=0.001
        sys.stdout.flush()
    
    walkers=survivors
    if len(survivors)<initialwalkernum:
        walkers = revive(survivors)
    survivors=list() #Note: this comes after the revival since old/new keep revivals separate
    
    thewalkingdead=0
    currAtZero = 0 
    for v in walkers: # move walkers, one walker at a time
        if (hammingWeight(v)==0):
            currAtZero+=1
            
        probedge= (1-s)*deltaT
        probdead= potential(v)*deltaT*s
        probstay = 1- probdead - probedge
        
        randprob=random.random()        
        if(randprob<=probedge): # walk to a neighbor
            neighbors=incidencedict[v]
            randvertex=random.choice(neighbors)
            survivors.append(randvertex)
        elif(randprob <= probedge+probdead): # walk to the cemetery
            thewalkingdead=thewalkingdead+1
        else: # stay
            survivors.append(v)
            
    if len(survivors) == 0:
        print("All walkers died")
        break
        
    writer.writerow((histogram(survivors)[0]))
    distfile.write(" %s \n" % (histogram(survivors)[0]))
    #if((t % int(1/deltaT))==0):
    #distfile.write(" %s %s %s " % (t,walkerlistnew[layerdict[0][0]], walkerlistnew[layerdict[1][0]]) )
    #distfile.write(" %s %s %s \n" % (walkerlistnew[layerdict[2][0]], walkerlistnew[layerdict[3][0]], thewalkingdead) )
    #distfile.write(" %s \n" % (walkerlistnew) )

time2=time.time()-time1
flag = False
for v in survivors:
    if potential(v)==0:
        flag=True
        break

if flag:
    print("\n Solution Found")
else:
    print("\n Failed to find solution")
#infofile.write(" %s \n" % (hammingdist) )
#infofile.write(" %s \n" % (vertexlist) )
#infofile.write(" %s %s \n" % (walkerlistold, thewalkingdead) )
#infofile.write(" %s %s %s %s \n" % (walkerlistnew, probedge, probdead, probstay) )
#infofile.write(" %s  \n" % (time2) )
