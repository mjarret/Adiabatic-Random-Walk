from sys import argv
import math
import random
import scipy
import networkx as nx
import time
script, filename1, filename2 = argv
print("usage: python rwmc.py info_file distribution_file")
infofile = open(filename1, 'w') # writing file
infofile.truncate()
distfile = open(filename2, 'w')  # writing file
distfile.truncate()
##################### Parameters
T=10
deltaT=0.01
timesteps=int(T/deltaT) # number of steps
n=10 # dimension of hypercube
if (deltaT > (1/float(3*n))): print("Warning: time-step deltaT is very large") # NOTE this condition depends on the dynamics
vertexnumber=2**n # for hypercube
initialwalkernum=4*n # initial number of walkers should be logarithmic in the number of vertices
##################### Use NetworkX to generate a hypercube; generate various dictionaries with the output
hyper=nx.hypercube_graph(n) # generate hypercube using NetworkX
vertexlist=hyper.nodes()
edgelist=hyper.edges()
adjlist=hyper.adjacency_list()
incidencedict=dict() # adjacency dictionary
hammingdist=dict() # hamming weight dictionary
layerdict=dict() # dictionary of vertices organized according to hamming weights
#eg for n=3 layerdict[2]={(1,1,0),(1,0,1),(0,1,1)}
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
walkerlistold=list() # list of walker locations
walkerlistnew=list()
for j in range(initialwalkernum): # randomize initial distribution
	rando = random.randint(0,vertexnumber-1)
	walkerlistnew.append(vertexlist[rando])
thewalkingdead=0 # initialize number of walkers that died in the last step
time1=time.time() # for timing of MCMC
for t in range(timesteps): # random walk loop	
	s=t/float(timesteps) # adiabatic scheduling parameter
	walkernum=0 # reset walkernum, to be recalculated and then called in the next two loops
	walkerlistold=list() # reset old list, to be updated in the following lines	
	for u in walkerlistnew:
		walkerlistold.append(u)
		walkernum=walkernum+1
	for count in range(thewalkingdead): # revival step # Note: range(0)=[]
		rando=random.randint(0,walkernum-1)
		w=walkerlistnew[rando] # pick a walker location randomly selected from the surviving walker distribution
		walkerlistold.append(w) # revive at w
	walkerlistnew=list() #Note: this comes after the revival since old/new keep revivals separate
	thewalkingdead=0 # reset the number of walkers sent to the cemetery
	for v in walkerlistold: # move walkers, one walker at a time
		probedge=deltaT*(1-s)*n
		probdead= deltaT*s*hammingdist[v]
		probstay=1-probedge-probdead
		randprob=random.random()
		if(randprob<=probedge): # walk to a neighbor
			randindex=random.randint(0,n-1) 
			neighbors=incidencedict[v]
			randvertex=neighbors[randindex]
			walkerlistnew.append(randvertex)
		elif((probedge< randprob) and (randprob <= probedge+probdead)): # walk to the cemetery
			thewalkingdead=thewalkingdead+1
		else: # stay
			walkerlistnew.append(v)
	#if((t % int(1/deltaT))==0):
	#distfile.write(" %s %s %s " % (t,walkerlistnew[layerdict[0][0]], walkerlistnew[layerdict[1][0]]) )
	#distfile.write(" %s %s %s \n" % (walkerlistnew[layerdict[2][0]], walkerlistnew[layerdict[3][0]], thewalkingdead) )
	#distfile.write(" %s \n" % (walkerlistnew) )
time2=time.time()-time1
#infofile.write(" %s \n" % (hammingdist) )
#infofile.write(" %s \n" % (vertexlist) )
#infofile.write(" %s %s \n" % (walkerlistold, thewalkingdead) )
#infofile.write(" %s %s %s %s \n" % (walkerlistnew, probedge, probdead, probstay) )
infofile.write(" %s  \n" % (time2) )
