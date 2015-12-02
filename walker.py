import random
import gmpy
import copy

def flipBit(number,  bit):
    return number^(2**bit)

def randomBitFlip(number, max):
    return flipBit(number, random.randint(0, max-1))
    
class Walker:
    def __init__(self,vertex=None, dimension=None):
       self.vertex=vertex
       self.dimension=dimension

    def getVertex(self):
        return self.vertex
        
    def walk(self):
        self.vertex=randomBitFlip(self.vertex, self.dimension)
        return self.vertex
        
    def potential(self):
        return gmpy.popcount(self.vertex)
        
    def spawn(self):
        return copy.copy(self)
        
