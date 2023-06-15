import time
import numpy as np
import matplotlib.pyplot as plt


# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getPairwiseSeparations(ri, rj):
    M = ri.shape[0]
    N = rj.shape[0]

    # positions ri = (x,y)
    rix = ri[:, 0].reshape((M, 1))
    riy = ri[:, 1].reshape((M, 1))

    # other set of points positions rj = (x,y)
    rjx = rj[:, 0].reshape((N, 1))
    rjy = rj[:, 1].reshape((N, 1))

    # matrices that store all pairwise particle separations: r_i - r_j
    dx = rix - rjx.T
    dy = riy - rjy.T

    return dx, dy


# set position of particles in a cube
def getPositionAndColorInGrid(N, nParticlesX, nParticlesY, distance):
    # top-left position of particle cube
    xPos = (nParticlesX * distance) / 2  # center the cube on x-axis

    pos = np.zeros((N, 2))

    x = np.tile(np.arange(nParticlesX), nParticlesY)
    y = np.repeat(np.arange(nParticlesY), nParticlesX)

    pos[:, 0] = x * distance - xPos
    pos[:, 1] = y * distance

    colors = np.full(N, "yellow")

    for i in range(N):
        if pos[i, 0] < 0:
            pos[i, 0] -= 3
            if pos[i, 1] > 10:
                pos[i, 1] += 10
                colors[i] = "red"
            else:
                colors[i] = "blue"

        else:
            pos[i, 0] += 3
            if pos[i, 1] > 10:
                pos[i, 1] += 10
                colors[i] = "yellow"
            else:
                colors[i] = "green"

    return pos, colors


# kernel function / W function
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def W(x, y, h):
    # kernel radius, include all surrounding particles of one particle in a grid
    r = np.sqrt(x**2 + y**2)
    return (1.0 / (h * np.sqrt(np.pi))) ** 3 * np.exp(-(r**2) / h**2)


# gradient of the kernel function
def gradW(x, y, h):
    r = np.sqrt(x**2 + y**2)

    n = -2 * np.exp(-(r**2) / h**2) / h**5 / (np.pi) ** (3 / 2)

    wx = n * x
    wy = n * y

    return wx, wy


# calculate density of all particles
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getDensity(mass, N, pos, h):
    dx, dy = getPairwiseSeparations(pos, pos)
    density = np.sum(mass * W(dx, dy, h), 1).reshape((N, 1))
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

    pressure = B * (((density / densityReference) ** m) - 1) + pressureReference
    pressure = np.maximum(pressure, 0)
    return pressure  # Pa


# calculate acceleration of all particles
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getAcceleration(N, mass, pos, vel, densityReference, h, nu):
    density = getDensity(mass, N, pos, h)
    pressure = getPressure(density, densityReference)

    gravity = (0, -1)

    dx, dy = getPairwiseSeparations(pos, pos)

    # Calculate the gradient of the kernel function
    dWx, dWy = np.array(gradW(dx, dy, h))

    # Calculate the pressure gradient
    ax = np.sum(
        -mass * (pressure / (density.T**2) + pressure.T / (density.T**2)) * dWx, 1
    ).reshape((N, 1))
    ay = np.sum(
        -mass * (pressure / (density.T**2) + pressure.T / (density.T**2)) * dWy, 1
    ).reshape((N, 1))

    acc = np.hstack((ax, ay))

    # gravity
    acc += gravity

    # viscosity
    acc -= nu * vel
    return acc


def main():
    # Simulation parameters

    # size of particle cube
    nParticlesX = 40  # number of particles in x direction
    nParticlesY = 40  # number of particles in y direction

    # distance between Particles at start
    distance = 0.5  # m

    # kernel radius
    h = np.sqrt(2) * distance  # m
    # number of particles
    N = nParticlesX * nParticlesY

    densityReference = 1000  # kg/m^2

    nu = 0.5

    # mass particle
    mass = densityReference * distance**2  # kg

    # initial positions of particles
    pos, colors = getPositionAndColorInGrid(N, nParticlesX, nParticlesY, distance)  # m
    # initial acceleration of particles
    acc = np.zeros((N, 2))  # m/s**2
    # initial velocity of particles
    vel = np.zeros((N, 2))  # m/s

    # show window
    fig, ax = plt.subplots()

    maxT = 10000
    dt = 0.1

    plt.show(block=False)

    # main loop
    for t in range(maxT):
        # (1/2) kick
        vel += acc * dt / 2

        # drift
        pos += vel * dt

        # collisions
        boundaryLeft = -20
        boundaryRight = 20
        boundaryTop = 100
        boundaryBottom = -3

        boundDamping = -0.8

        vel[(pos[:, 0] > boundaryRight), 0] *= boundDamping
        pos[(pos[:, 0] > boundaryRight), 0] = boundaryRight

        vel[(pos[:, 0] < boundaryLeft), 0] *= boundDamping
        pos[(pos[:, 0] < boundaryLeft), 0] = boundaryLeft

        vel[(pos[:, 1] > boundaryTop), 1] *= boundDamping
        pos[(pos[:, 1] > boundaryTop), 1] = boundaryTop

        vel[(pos[:, 1] < boundaryBottom), 1] *= boundDamping
        pos[(pos[:, 1] < boundaryBottom), 1] = boundaryBottom

        # update accelerations
        acc = getAcceleration(N, mass, pos, vel, densityReference, h, nu)

        # (1/2) kick
        vel += acc * dt / 2

        # update window
        plt.sca(ax)
        plt.cla()
        plt.scatter(
            pos[:, 0],
            pos[:, 1],
            s=20,
            alpha=0.5,
            c=colors,
        )
        ax.set(xlim=(-20, 20), ylim=(-6, 20))
        plt.pause(0.0001)

    return 0


if __name__ == "__main__":
    main()
