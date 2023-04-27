import numpy as np
from math import pi
import random

#number of particles
N = 70
#mass of particles
m = 1
#kernel radius
h = 100
#density of particles
d = np.zeros(N)
#positions of particles
pos = np.zeros((N, 2))
#size of the screen
size = 1000

p = np.zeros(N)

#setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random()*size;
        pos[i][1] = random.random()*size;

#W function
def W(r, h):
    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))

#calculate density of all particles
def calculateDensity():
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.linalg.norm(pos[i]-pos[j])
                d[i] += m * W(r,h)

#calculate pressure with density of all particles
def calculatePressure():
    for i in range (N):
        for j in range (N):
            if (i!=j):
                p[i] += 0.3214 * ((d[i]) ** 7 - 1) + 1

def main():
    setRandomPosition()
    calculateDensity()
    calculatePressure()
    print(p)

main()
print()
print(d)