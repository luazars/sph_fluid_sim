import numpy as np
import pygame
from math import pi
import random

N = 1000
m = 1
size = 1000
h = 100

pos = np.zeros((N, 2))
for i in range(N):
    pos[i][0] = random.random()*size;
    pos[i][1] = random.random()*size;

p = np.zeros(N)

def W(r, h):
    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))


def P():
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.linalg.norm(pos[i]-pos[j])
                p[i] += m * W(r,h)
    

pygame.init()


screen = pygame.display.set_mode([size,size])

done = False
clock = pygame.time.Clock()

P()

while not done:
    clock.tick(1) 
    for event in pygame.event.get():  
        if event.type == pygame.QUIT:
            done = True

    screen.fill("white")
    for i in range(N):
        pygame.draw.circle(screen, [ (p[i]/max(p))*255,0,0 ], [pos[i][0], pos[i][1]], 5)

    pygame.display.flip()

pygame.quit()

