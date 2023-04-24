import numpy as np
import pygame
from math import pi
import random

#number of particles
N = 1000
#mas of particles
m = 1
#kernel radius
h = 100
#density of particles
d = np.zeros(N)
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
def W(r, h):
    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))

#calculate density of all particles
def calculateDensity():
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.linalg.norm(pos[i]-pos[j])
                d[i] += m * W(r,h)

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

        #drawing each particle
        for i in range(N):
            pygame.draw.circle(screen, [ (d[i]/max(d))*255,0,0 ], [pos[i][0], pos[i][1]], 5)

        #update screen
        pygame.display.flip()

    pygame.quit()


def main():
    setRandomPosition()
    calculateDensity()
    showWindow()


if __name__ == "__main__":
    main()