import numpy as np
import pygame
from math import pi
import random

#width and height of the screen
size = 1000 #mm

#size of particle cube
xSize = 50 #mm
ySize = 50 #mm

#distance between Particles at start
distance = 2 #mm

#number of particles
N = xSize*ySize

#density of particles
density = np.zeros(N) #u/mm**2
#pressure of particles
pressure = np.zeros(N) #MPa
#positions of particles
pos = np.zeros((N, 2)) #mm

#setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random()*size;
        pos[i][1] = random.random()*size;

#set position of particles in a cube
def setPositionInGrid():
    #top-left position of particle cube
    xPos = 1/2 * (size-xSize*distance) #center the cube on x-axis 
    yPos = 100

    for i in range(N):
        x = i%xSize
        y =  int(i/N*ySize)
        pos[i][0] = x*distance+xPos
        pos[i][1] = y*distance+yPos

#kernel function / W function
def W(r):
    #kernel radius, include all surrounding particles of one particle in a grid
    h = 2**(1/2)*distance #mm

    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))

#calculate density of all particles
def calculateDensity():

    #mass of water particles
    mass =  18.01528 #u

    for i in range(N):
        for j in range(N):
            r = np.linalg.norm(pos[i]-pos[j])
            density[i] += mass * W(r)
        

#calculate pressure with density of all particles with tait equation
def calculatePressure():
    #?
    B = 321.4 #MPa
    #Adiabatenexponent
    m = 7 

    pressureReference = 2060.0 #MPa
    densityReference = 6.004  #komische einheit
    
    for i in range (N):
        pressure[i] = pressureReference + B * ((density[i] / densityReference) ** m - 1) #MPa
        print(pressure[i])

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
            pygame.draw.circle(screen, [(density[i]-min(density)) /(max(density)-min(density))*255,0,0 ], [pos[i][0], pos[i][1]], ((pressure[i]-min(pressure)) /(max(pressure)-min(pressure)))*5+4)

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