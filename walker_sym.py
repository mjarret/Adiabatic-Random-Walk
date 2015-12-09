import random
import gmpy
import copy

def flipBit(number,  bit):
    return number^(2**bit)

def randomBitFlip(number, max):
    return flipBit(number, random.randint(0, max-1))
    
class Walker:
    def __init__(self,vertex=None, dimension=None):
       self.hw = gmpy.popcount(vertex)
       self.dimension=dimension

    def getVertex(self):
        return self.vertex
        
    def walk(self):
        down_steps = self.hw
        up_steps = self.dimension-self.hw
        r = random.randrange(down_steps + up_steps)
        if(r < down_steps): self.hw -=1
        else: self.hw +=1      
        return self.hw
      
    def hammingWeight(self):
        return gmpy.popcount(self.vertex)
        
    def spawn(self):
        return copy.copy(self)
        
