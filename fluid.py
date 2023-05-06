import numpy as np
import pygame
from math import pi
import random

#size of the screen
size = 1000

#distance between Particles at start
distance = 15

#size of particle cube
xSize = 50
ySize = 10

#pos of particle cube
xPos = 1/2 * (size-xSize*distance) #center the cube on x-axis 
yPos = 100

#number of particles
N = xSize*ySize

#mass of particles
mass = 1

#kernel radius, include all surrounding particles of one particle in a grid
h = 2**(1/2)*distance

#?
B = 5.0
#Adiabatenexponent
m = 2.5

startingPressure = 0.0
startingDensity = 1.0

#density of particles
density = np.zeros(N)
#pressure of particles
pressure = np.zeros(N)
#positions of particles
pos = np.zeros((N, 2))

#setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random()*size;
        pos[i][1] = random.random()*size;

#set position of particles in a cube
def setPositionInGrid():
    for i in range(N):
        x = i%xSize
        y =  int(i/N*ySize)
        pos[i][0] = x*distance+xPos
        pos[i][1] = y*distance+yPos

#kernel function / W function
def W(r):
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
    for i in range (N):
        pressure[i] = startingPressure + B * ((density[i] / startingDensity) ** m - 1)   

def showWindow():
    #initialization of pygame 
    pygame.init()
    screen = pygame.display.set_mode([size,size])
    done = False
    clock = pygame.time.Clock()

    #drawing loop
    while not done:
        clock.tick(1)

        #checking if closing button is pressed
        for event in pygame.event.get():  
            if event.type == pygame.QUIT:
                done = True

        #background color
        screen.fill("white")

        #drawing each particle with color according to the density and radius according to pressure
        for i in range(N):
            pygame.draw.circle(screen, [ (density[i]-min(density)) /(max(density)-min(density))*255,0,0 ], [pos[i][0], pos[i][1]], ((pressure[i]-min(pressure)) /(max(pressure)-min(pressure)))*5+4)

        #update screen
        pygame.display.flip()

    pygame.quit()


def main():
    setPositionInGrid()
    calculateDensity()
    calculatePressure()
    showWindow()


if __name__ == "__main__":
    main()