import sys
#import math
import random
import numpy
import argparse
import walker
import csv
import gmpy


# Some parsers for command line options
parser = argparse.ArgumentParser(description='Walk with killing.')
parser.add_argument('infofile',  metavar='infofile',  nargs='?',  help='output file 1',  default='info.csv')
parser.add_argument('distfile',  metavar='distfile',  nargs='?',  help='output file 2',  default='dist')
parser.add_argument('-T', metavar='T',  nargs='?',  type=int,  help='total runtime', default=100)
parser.add_argument('-s', metavar='s',  nargs='?', type=float,  help='discretization size', default=.3)
parser.add_argument('-w', metavar='w',  nargs='?', type=int,  help='total number of walkers', default=100)
parser.add_argument('-d', metavar='d',  nargs='?', type=int,  help='dimension of hypercube graph', default=10)
parser.add_argument('--spike_size', metavar="spike_size",  nargs='?', type=float, help='introduce a spike of height spike_size',  default=0)
parser.add_argument('--spike_loc', metavar="spike_loc",  nargs='?', type=int, help='place spike at hamming weight spike_loc',  default=0)
args=parser.parse_args()
#print("usage: python adiabaticv2.py info_file distribution_file")
infofile = open(args.infofile, 'w') # writing file
infofile.truncate()
distfile = open(args.distfile, 'w')  # writing file
distfile.truncate()
n=args.d # dimension of hypercube
initialwalkernum=args.w # initial number of walkers should be logarithmic in the number of verticesu
deltaT=args.s
vertexnumber=int(2**n) # for hypercube
writer = csv.writer(infofile)

spike_size = args.spike_size
spike_loc = args.spike_loc

def main():
    #script, filename1, filename2 = argv
    T = args.T
 #   results = list()
    for time_counter in range(10):
        w=args.w
        while(w < 10*(n**2)):
            trials = 5
            tmp = list()
            for i in range(trials):
                tmp.append(test(T, w))
            tmp[:0] = (n, T, w, spike_loc, spike_size)
            writer.writerow(tmp)
            print(tmp)
            w = 2*w
        T = 2*T
        
    
    
##################### Use NetworkX to generate a hypercube; generate various dictionaries with the output
       
#def revive( walkerList ):
#    tmpList = list()
#    tmpList.extend(walkerList)
#    while len(tmpList) < initialwalkernum:
#        w = random.choice(walkerList)                           # choose a random surviving member
#        tmpList.append(w.spawn())                            # revive one
#    return tmpList
    
binList = list(range(n+1))
def histogram( walkerList ):
    histoList = list()
    for w in walkerList:
        histoList.append(w.hammingWeight())
    return numpy.histogram(histoList, bins=binList, density=False)

def adiabaticWalk(total_walkers, total_vertices, time):
    timesteps=int(time/deltaT)
#    percentDone = 0.000
    walkers = list()
    for j in range(total_walkers):           # randomize initial distribution
        r = random.randrange(total_vertices)
        walkers.append(walker.Walker(r, total_vertices)) #add a new walker

    for t in range(timesteps): # random walk loop
        s=t/float(timesteps) # adiabatic scheduling parameter
    
    #    writer.writerow((histogram(walkers)[0]))
        distfile.write(" %s \n" % (histogram(walkers)[0]))
    
#        if s > percentDone:
#            sys.stdout.write("\r %s %% Complete" % (percentDone*100))
#            percentDone+=0.001
#            sys.stdout.flush()

        walkers = diffuse(walkers, s)
    return walkers

def diffuse(prior, s):   
    survivors = list()          
    for w in prior:
            
        probedge= (1-s)*deltaT
        probdead= potential(w.vertex)*deltaT*s
        probstay = 1- probdead - probedge
        
        randprob=random.random()        
        
        if(randprob<=probedge): 
            w.walk()
            survivors.append(w)
        elif(randprob <= probedge+probstay): 
            survivors.append(w)
        else:
            survivors.append(random.choice(prior).spawn())
    return survivors 

def test(time,  num_walkers):    
    walkers = adiabaticWalk(num_walkers, n, time)
    hit = 0
    for w in walkers:
        if w.vertex == 0:
            hit += 1
    return hit/float(num_walkers)

def hammingWeight(vertex):
        return gmpy.popcount(vertex)
        
def potential(vertex):
    h=hammingWeight(vertex)
    if (h == spike_loc): return spike_size/float(n)
    return h/float(n)

if __name__ == "__main__":
    sys.exit(main())
