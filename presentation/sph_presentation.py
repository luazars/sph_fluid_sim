import time
import numpy as np
import matplotlib.pyplot as plt


def getPositionInGrid(N, nParticlesX, nParticlesY, distance):
    pos = np.zeros((N, 2))

    xPos = (nParticlesX * distance) / 2  # center the cube on x-axis

    x = np.tile(np.arange(nParticlesX), nParticlesY)
    y = np.repeat(np.arange(nParticlesY), nParticlesX)

    pos[:, 0] = x * distance - xPos
    pos[:, 1] = y * distance + 2

    return pos


def W(r, h):
    return (1.0 / (h * np.sqrt(np.pi))) ** 3 * np.exp(-(r**2) / h**2)


def gradW(x, y, r, h):
    n = -2 * np.exp(-(r**2) / h**2) / h**5 / (np.pi) ** (3 / 2)
    wx = n * x
    wy = n * y
    return wx, wy


def getDensity(mass, N, pos, h):
    density = np.zeros(N)

    for i in range(N):
        for j in range(N):
            r = np.linalg.norm(pos[i] - pos[j])
            density[i] += mass * W(r, h)

    return density


def getPressure(density, densityReference, pressureReference):
    c = 5  # Speed of sound
    m = 7  # Adiabatenexponent?
    B = (densityReference * c**2) / m

    pressure = B * (((density / densityReference) ** m) - 1) + pressureReference
    pressure = np.maximum(pressure, 0)

    return pressure


def getAcceleration(N, mass, pos, vel, densityReference, pressureReference, h, nu, g):
    density = getDensity(mass, N, pos, h)
    pressure = getPressure(density, densityReference, pressureReference)
    acc = np.zeros((N, 2))

    for i in range(N):
        for j in range(N):
            if i != j:
                r = pos[i] - pos[j]
                dist = np.linalg.norm(r)

                grad = np.array(gradW(r[0], r[1], dist, h))

                pressure_grad = (
                    -mass
                    * (
                        pressure[i] / (density[i] ** 2)
                        + pressure[j] / (density[j] ** 2)
                    )
                    * grad
                )

                acc[i] += pressure_grad

        # gravity
        acc[i] += g

        # viscosity
        acc[i] -= nu * vel[i]

    return acc


def main():
    # * Simulation parameters

    # size of particle cube
    nParticlesX = 10  # number of particles in x direction
    nParticlesY = 10  # number of particles in y direction

    # number of particles
    N = nParticlesX * nParticlesY

    # distance between Particles at start
    distance = 0.5  # m

    # kernel radius
    h = np.sqrt(2) * distance  # m

    densityReference = 700  # kg/m^2
    pressureReference = 0

    # mass particle
    mass = densityReference * distance**2  # kg

    # weight force
    g = (0, -2)

    # viscosity
    nu = 0.2

    # initial positions of particles
    pos = getPositionInGrid(N, nParticlesX, nParticlesY, distance)  # m
    # initial velocity of particles
    vel = np.zeros((N, 2))  # m/s
    # initial acceleration of particles
    acc = getAcceleration(
        N, mass, pos, vel, densityReference, pressureReference, h, nu, g
    )  # m/s**2

    maxT = 200
    dt = 0.06

    # show window
    fig, ax = plt.subplots()
    plt.show(block=False)

    # * main loop
    for t in range(maxT):
        # (1/2) kick
        vel += acc * dt / 2

        # drift
        pos += vel * dt

        # boundary positions
        boundaryLeft = -4
        boundaryRight = 4
        boundaryTop = 1000
        boundaryBottom = 0

        boundaryDamping = -0.5

        for i in range(N):
            # right boundary
            if pos[i][0] > boundaryRight:
                vel[i][0] *= boundaryDamping
                pos[i][0] = boundaryRight
            # left boundary
            if pos[i][0] < boundaryLeft:
                vel[i][0] *= boundaryDamping
                pos[i][0] = boundaryLeft

            # top boundary
            if pos[i][1] > boundaryTop:
                vel[i][1] *= boundaryDamping
                pos[i][1] = boundaryTop
            # bottom boundary
            if pos[i][1] < boundaryBottom:
                vel[i][1] *= boundaryDamping
                pos[i][1] = boundaryBottom

        # update accelerations
        acc = getAcceleration(
            N, mass, pos, vel, densityReference, pressureReference, h, nu, g
        )

        # (1/2) kick
        vel += acc * dt / 2

        # update window
        plt.sca(ax)
        plt.cla()
        plt.scatter(
            pos[:, 0],
            pos[:, 1],
            s=40,
            alpha=0.5,
        )
        ax.set(xlim=(-4, 4), ylim=(-0, 8))
        plt.gca().set_aspect("equal")
        plt.pause(0.001)

    return 0


if __name__ == "__main__":
    main()
