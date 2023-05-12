import numpy as np
import pygame
from math import pi
import random

# width and height of the screen
size = 1  # m

# size of particle cube
nParticlesX = 40  # number of particles in x direction
nParticlesY = 40  # number of particles in y direction

# distance between Particles at start
distance = 0.05  # m

# number of particles
N = nParticlesX * nParticlesY


# speed of sound
c = 15  # m/s
# Adiabatenexponent
m = 7

densityReference = 1000  # kg/m^2

# mass particle
mass = densityReference * distance**2  # kg

# density of particles
density = np.zeros(N)  # kg/m^2
# pressure of particles
pressure = np.zeros(N)  # Pa
# positions of particles
pos = np.zeros((N, 2))  # m


# setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random() * size
        pos[i][1] = random.random() * size


# set position of particles in a cube
def setPositionInGrid():
    # top-left position of particle cube
    xPos = 1 / 2 * (size - nParticlesX * distance)  # center the cube on x-axis
    yPos = 0.1

    for i in range(N):
        x = i % nParticlesX
        y = int(i / N * nParticlesY)
        pos[i][0] = x * distance + xPos
        pos[i][1] = y * distance + yPos


# kernel function / W function
def W(r):
    # kernel radius, include all surrounding particles of one particle in a grid
    h = 2 ** (1 / 2) * distance  # m

    return (1 / ((h**3) * (np.pi ** (2 / 3)))) * np.exp((-(np.abs(r) ** 2) / h**2))


# calculate density of all particles
def calculateDensity():
    for i in range(N):
        for j in range(N):
            r = np.linalg.norm(pos[i] - pos[j])
            density[i] += mass * W(r)


# calculate pressure with density of all particles with tait equation / Zustandsgleichung
def calculatePressure():
    for i in range(N):
        pressure[i] = (
            (densityReference * c**2) / m * ((density[i] / densityReference) ** m - 1)
        )  # Pa


# matplotlib / pyplot /scatterplot
def showWindow():
    # initialization of pygame
    pygame.init()
    screen = pygame.display.set_mode([size * 1000, size * 1000])
    done = False
    clock = pygame.time.Clock()

    # drawing loop
    while not done:
        clock.tick(1)

        # checking if closing button is pressed
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                done = True

        # background color
        screen.fill("white")

        # drawing each particle with color according to the density and radius according to pressure
        for i in range(N):
            pygame.draw.circle(
                screen,
                [
                    (density[i] - min(density)) / (max(density) - min(density)) * 255,
                    0,
                    0,
                ],
                [pos[i][0] * 1000, pos[i][1] * 1000],
                ((pressure[i] - min(pressure)) / (max(pressure) - min(pressure))) * 10
                + 1,
            )

        # update screen
        pygame.display.flip()

    pygame.quit()


def main():
    setPositionInGrid()
    calculateDensity()
    calculatePressure()
    showWindow()


if __name__ == "__main__":
    main()
