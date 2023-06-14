import time
import numpy as np
import matplotlib.pyplot as plt


# set position of particles in a cube
def getPositionInGrid(N, nParticlesX, nParticlesY, distance):
    # top-left position of particle cube
    xPos = (nParticlesX * distance) / 2  # center the cube on x-axis

    pos = np.zeros((N, 2))
    pos = np.zeros((N, 2))
    x = np.tile(np.arange(nParticlesX), nParticlesY)
    y = np.repeat(np.arange(nParticlesY), nParticlesX)

    pos[:, 0] = x * distance - xPos
    pos[:, 1] = y * distance
    return pos


# kernel function / W function
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def W(r, h):
    # kernel radius, include all surrounding particles of one particle in a grid
    return (1.0 / (h * np.sqrt(np.pi))) ** 3 * np.exp(-(r**2) / h**2)


# gradient of the kernel function
def gradW(x, y, r, h):
    n = -2 * np.exp(-(r**2) / h**2) / h**5 / (np.pi) ** (3 / 2)
    wx = n * x
    wy = n * y
    return wx, wy


# calculate density of all particles
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getDensity(i, mass, N, pos, h):
    density = 0
    for j in range(N):
        r = np.linalg.norm(pos[i] - pos[j])
        density += mass * W(r, h)
    return density


# calculate pressure with density of all particles with tait equation
# https://de.wikipedia.org/wiki/Taitsche_Gleichung
def getPressure(density, densityReference):
    # speed of sound
    c = 5  # m/s
    # Adiabatenexponent
    m = 7
    pressureReference = 0

    B = (densityReference * c**2) / m
    return max(
        0, B * (((density / densityReference) ** m) - 1) + pressureReference
    )  # Pa


# calculate acceleration of all particles
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getAcceleration(N, mass, pos, density, pressure, h):
    # reset acceleration for particles
    acc = np.zeros((N, 2))
    for i in range(N):
        for j in range(N):
            if i != j:
                r = pos[i] - pos[j]
                dist = np.linalg.norm(r)

                # Calculate the gradient of the kernel function
                grad = np.array(gradW(r[0], r[1], dist, h))

                # Calculate the pressure gradient
                pressure_grad = (
                    -mass
                    * (
                        pressure[i] / (density[i] ** 2)
                        + pressure[j] / (density[j] ** 2)
                    )
                    * grad
                )

                acc[i] += pressure_grad

    acc[:, 1] += -2
    return acc


def timeStep(
    N, mass, dt, acc, vel, pos, density, densityReference, pressure, h, ax, fig
):
    # (1/2) kick
    vel += acc * dt / 2

    # drift
    pos += vel * dt

    for i in range(N):
        density[i] = getDensity(i, mass, N, pos, h)
        pressure[i] = getPressure(density[i], densityReference)

    # update accelerations
    acc = getAcceleration(N, mass, pos, density, pressure, h)

    # (1/2) kick
    vel += acc * dt / 2

    boundaryLeft = -5
    boundaryRight = 5
    boundaryTop = 10
    boundaryBottom = -5

    boundDamping = -0.5

    vel[(pos[:, 0] > boundaryRight), 0] *= boundDamping
    pos[(pos[:, 0] > boundaryRight), 0] = boundaryRight

    vel[(pos[:, 0] < boundaryLeft), 0] *= boundDamping
    pos[(pos[:, 0] < boundaryLeft), 0] = boundaryLeft

    vel[(pos[:, 1] > boundaryTop), 1] *= boundDamping
    pos[(pos[:, 1] > boundaryTop), 1] = boundaryTop

    vel[(pos[:, 1] < boundaryBottom), 1] *= boundDamping
    pos[(pos[:, 1] < boundaryBottom), 1] = boundaryBottom


def main():
    # Simulation parameters

    # size of particle cube
    nParticlesX = 10  # number of particles in x direction
    nParticlesY = 20  # number of particles in y direction

    # distance between Particles at start
    distance = 0.3  # m

    # kernel radius
    h = np.sqrt(2) * distance  # m
    # number of particles
    N = nParticlesX * nParticlesY

    densityReference = 1000  # kg/m^2

    # mass particle
    mass = densityReference * distance**2  # kg

    # initial positions of particles
    pos = getPositionInGrid(N, nParticlesX, nParticlesY, distance)  # m

    # density of particles
    density = np.zeros(N)  # kg/m^2
    # pressure of particles
    pressure = np.zeros(N)  # Pa

    # get initial density and pressure
    for i in range(N):
        density[i] = getDensity(i, mass, N, pos, h)
        pressure[i] = getPressure(density[i], densityReference)

    # initial acceleration of particles
    acc = getAcceleration(N, mass, pos, density, pressure, h)  # m/s**2

    # initial velocity of particles
    vel = np.zeros((N, 2))  # m/s

    # show window
    fig, ax = plt.subplots()
    plt.ion()

    maxT = 10000
    dt = 0.05

    # main loop
    for t in range(maxT):
        timeStep(
            N, mass, dt, acc, vel, pos, density, densityReference, pressure, h, ax, fig
        )
        plt.sca(ax)
        plt.cla()
        plt.scatter(pos[:, 0], pos[:, 1], s=10, c=density, alpha=0.5)
        ax.set(xlim=(-10, 10), ylim=(-10, 10))

        plt.pause(0.001)

    plt.show(block=False)
    plt.close("all")
    return 0


if __name__ == "__main__":
    main()
