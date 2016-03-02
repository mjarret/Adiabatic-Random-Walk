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
parser.add_argument('--post_select', metavar="post_select",  nargs='?',  type=bool,  help='condition on distribution after',  default = False)
args=parser.parse_args()
#print("usage: python adiabaticv2.py info_file distribution_file")

spike_size = args.spike_size
spike_loc = args.spike_loc
n=args.d # dimension of hypercube
curr_min = 100000
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

dim = args.w
post_select = args.post_select

tot_walkers = args.w

currfile = open('out/curr.csv', 'a+')
curr = csv.writer(currfile)

initialwalkernum=args.w # initial number of walkers should be logarithmic in the number of verticesu
deltaT=args.s
vertexnumber=int(2**n) # for hypercube

writer = csv.writer(infofile)
dist_writer = csv.writer(distfile)

spike_size = args.spike_size
spike_loc = args.spike_loc

stephen_b = 2/float(numpy.tan(2*numpy.arccos(1 - 1/float(n))))
print(stephen_b)

def main():
    T = args.T
    
    # These loops aren't typically used
    for time_counter in range(time_range):
        w=args.w
        while(w < 10*(n**2)):
            tmp = list()
            tmp = [test(T, w) for i in range(trials)]
            curr.writerow((n, T, w, spike_loc, spike_size,numpy.mean(tmp, dtype=numpy.float64)))
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
            output = [s]
            ll = histogram(walkers)[0].tolist()
            output.extend(ll)
            dist_writer.writerow(output)
            s_old += .001
        walkers = diffuse(walkers, s)
    
    # sawtooth...
#    q=0.01
#    s_old = 0.000
#    for m in my_range(q, 1, q):
#        for t in my_range((m-20*q)*timesteps, timesteps, 1): # random walk loop
#            s=t/float(timesteps) # adiabatic scheduling parameter
#            if(s < 0): break
#            if(s > 1): break
#            if(s - s_old > .001): 
#                output = [s]
#                ll = histogram(walkers)[0].tolist()
#                output.extend(ll)
#                dist_writer.writerow(output)
#                s_old += .001
#            walkers = diffuse(walkers, s)
    return walkers
    

def diffuse(prior, s):
    s = 1/float(2)   
    survivors = list() 
    av = 0.0
    n = len(prior)
    av = numpy.mean([potential(w.hw) for w in prior])
    diff = -(tot_walkers - n)/float(n)
    
    for w in prior:
        pot = potential(w.hw) - av + diff
        probedge = (1-s)*deltaT
        probdead = s*abs(pot)*deltaT
        probstay = 1- probdead - probedge

        randprob=random.random()        

        if(randprob<=probedge): 
            w.walk()
            survivors.append(w)
        elif(randprob <= probedge+probstay): 
            survivors.append(w)
        elif(pot < 0):
            survivors.append(w)
            survivors.append(w.spawn())
    return survivors
    
def my_range(start, end, step):
    while start <= end:
        yield start
        start += step    

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
    # Birth/Death calculates probabilities based on 2W - 1 for W between 0 and 1. Need to return 2W instead of W.
    if (hamming_weight == spike_loc): return 2*spike_size/float(n)
    h= hamming_weight
    pot = stephen_b*2*h/float(n)
    #pot = h/float(n)
    return pot


if __name__ == "__main__":
    sys.exit(main())
