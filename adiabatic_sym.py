import sys
#import math
import random
import numpy
import argparse
import walker_sym
import csv
import gmpy
import os


# Some parsers for command line options
parser = argparse.ArgumentParser(description='Walk with killing.')
parser.add_argument('infofile',  metavar='infofile',  nargs='?',  help='output file 1',  default='info.csv')
parser.add_argument('distfile',  metavar='distfile',  nargs='?',  help='output file 2',  default='dist.csv')
parser.add_argument('-T', metavar='T',  nargs='?',  type=int,  help='total runtime', default=100)
parser.add_argument('-s', metavar='s',  nargs='?', type=float,  help='discretization size', default=.3)
parser.add_argument('-w', metavar='w',  nargs='?', type=int,  help='total number of walkers', default=100)
parser.add_argument('-d', metavar='d',  nargs='?', type=int,  help='dimension of hypercube graph', default=10)
parser.add_argument('--spike_size', metavar="spike_size",  nargs='?', type=float, help='introduce a spike of height spike_size',  default=0)
parser.add_argument('--spike_loc', metavar="spike_loc",  nargs='?', type=int, help='place spike at hamming weight spike_loc',  default=0)
args=parser.parse_args()
#print("usage: python adiabaticv2.py info_file distribution_file")

spike_size = args.spike_size
spike_loc = args.spike_loc
n=args.d # dimension of hypercube

trials = 1
time_range = 1

dir = 'out/%s/%s/%s/' % (n, spike_loc, spike_size)
if not os.path.exists(dir):
    os.makedirs(dir)
outfile = '%s/%s' % (dir, args.infofile)    
distfile = '%s/%s' % (dir, args.distfile)
infofile = open(outfile, 'a+') # writing file
distfile = open(distfile, 'w')  # writing file
distfile.truncate()

initialwalkernum=args.w # initial number of walkers should be logarithmic in the number of verticesu
deltaT=args.s
vertexnumber=int(2**n) # for hypercube

writer = csv.writer(infofile)
dist_writer = csv.writer(distfile)

spike_size = args.spike_size
spike_loc = args.spike_loc
    
def main():
    T = args.T
    for time_counter in range(time_range):
        w=args.w
        while(w < 10*(n**2)):
            tmp = list()
            tmp = [test(T, w) for i in range(trials)]
            tmp[:0] = (n, T, w, spike_loc, spike_size)
            writer.writerow(tmp)
            w = 2*w
            break
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
    for w in walkerList: histoList.append(w.hw)
    return numpy.histogram(histoList, bins=binList, density=True)

def adiabaticWalk(total_walkers, total_vertices, time):
    timesteps=int(time/deltaT)
#    percentDone = 0.000
    walkers = list()
    walkers = [walker_sym.Walker(random.randrange(2**total_vertices), total_vertices) for j in range(total_walkers)]           # randomize initial distribution
    s_old = 0.000
    for t in range(timesteps): # random walk loop
        s=t/float(timesteps) # adiabatic scheduling parameter
        if(s - s_old > .001): 
            print(s)
            output = [s]
            ll = histogram(walkers)[0].tolist()
            output.extend(ll)
            dist_writer.writerow(output)
            s_old += .001
        before_diffusion = walkers
        walkers = diffuse(walkers, s)
#        if (len(walkers) < total_walkers/float(8)):
#            while(len(walkers) < total_walkers): walkers.append(random.choice(before_diffusion).spawn())
    return walkers
    

def diffuse(prior, s):   
    survivors = list()          
    
    for w in prior:
            
        #probedge= (1-s)*deltaT
        #probdead= potential(w.hw)*deltaT*s
        probedge = (1-s)*deltaT
        probdead = s*potential(w.hw)*deltaT
        probstay = 1- probdead - probedge
        
        randprob=random.random()        
        
        if(randprob<=probedge): 
            w.walk()
            survivors.append(w)
        elif(randprob <= probedge+probstay): 
            survivors.append(w)
        else:
            survivors.append(random.choice(prior).spawn())
#            pass

        #probedge= (1-s)*deltaT
        #probdead= potential(w.hw)*deltaT*s
#        isDead = False
#        for i in range(2):
#            probedge = (1-s)
#            probdead = s*potential(w.hw)
#            probstay = 1- probdead - probedge
#        
#            randprob=random.random()        
#        
#            if(randprob<=probedge): 
#                w.walk()
#            elif(randprob <= probedge+probstay): 
#                pass
#            else:isDead = True
#                
#        if (isDead):        
#            survivors.append(random.choice(prior).spawn())
#        else:
#            survivors.append(w)
    return survivors 
    

def test(time,  num_walkers):    
    walkers = adiabaticWalk(num_walkers, n, time)
    hit = 0
    for w in walkers:
        if w.hw == 0:
            hit += 1
    return hit/float(num_walkers)

def hammingWeight(vertex):
        return gmpy.popcount(vertex)
        
def potential(hamming_weight):
    h=hamming_weight
    if (h == spike_loc): return spike_size/float(n)
    return h/float(n)

if __name__ == "__main__":
    sys.exit(main())
