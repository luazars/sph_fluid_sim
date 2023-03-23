import numpy as np
import pygame
from math import pi

N = 1000
m = 1

pos = np.random.randn(N,3)*100
p = np.zeros(N)

def W(r, h):
    return (1 / (h * np.sqrt(2*np.pi))) * np.exp((-(r**2)/2*h**2))


def P():
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.linalg.norm(pos[i]-pos[j])
                p[i] += m*W(r,.01)
    

pygame.init()

size = 1000
screen = pygame.display.set_mode([size,size])

done = False
clock = pygame.time.Clock()

p = np.zeros(N)
pos = np.random.randn(N,3)*100
P()

while not done:
    clock.tick(1) 
    for event in pygame.event.get():  
        if event.type == pygame.QUIT:
            done = True

    screen.fill("white")

    for i in range(N):
        pygame.draw.circle(screen, [ (p[i]/max(p))*255,0,0 ], [pos[i][0]+500, pos[i][1]+500], 5)

    pygame.display.flip()

pygame.quit()

