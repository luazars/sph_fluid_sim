import numpy as np
import pygame
from math import pi
import random

#number of particles
N = 1000
#mass of particles
mass = 1
#density of particles
density = np.zeros(N)
#pressure of particles
pressure = np.zeros(N)
#positions of particles
pos = np.zeros((N, 2))

#size of the screen
size = 1000

#setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random()*size;
        pos[i][1] = random.random()*size;

#W function
def W(r):
    #kernel radius
    h = 100
    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))

#calculate density of all particles
def calculateDensity():
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.linalg.norm(pos[i]-pos[j])
                density[i] += mass * W(r)

#calculate pressure with density of all particles with tait equation
def calculatePressure():
    #random values for B and m
    B = 2.0
    m = 1.5
    startingPressure = 0.0
    startingDensity = 1.0
    for i in range (N):
        pressure[i] = startingPressure + B * ((density[i] / startingDensity) ** m - 1)
    print(pressure)    

def showWindow():
    #initialization of pygame 
    pygame.init()
    screen = pygame.display.set_mode([size,size])
    done = False
    clock = pygame.time.Clock()

    #drawing loop
    while not done:
        clock.tick(1)

        #checking if window is closed
        for event in pygame.event.get():  
            if event.type == pygame.QUIT:
                done = True

        #background color
        screen.fill("white")

        #drawing each particle with color according to the density and radius according to pressure
        for i in range(N):
            pygame.draw.circle(screen, [ (density[i]-min(density)) /(max(density)-min(density))*255,0,0 ], [pos[i][0], pos[i][1]], ((pressure[i]-min(pressure)) /(max(pressure)-min(pressure)))*10+4)

        #update screen
        pygame.display.flip()

    pygame.quit()


def main():
    setRandomPosition()
    calculateDensity()
    calculatePressure()
    showWindow()


if __name__ == "__main__":
    main()