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


def getPressure(density, densityReference):
    c = 5
    m = 7
    pressureReference = 0
    B = (densityReference * c**2) / m

    pressure = B * (((density / densityReference) ** m) - 1) + pressureReference
    pressure = np.maximum(pressure, 0)
    return pressure


def getAcceleration(N, mass, pos, vel, densityReference, h, nu):
    density = getDensity(mass, N, pos, h)
    pressure = getPressure(density, densityReference)
    acc = np.zeros((N, 2))

    gravity = (0, -2)

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

        acc[i] += gravity
        acc[i] -= nu * vel[i]

    return acc


def main():
    # Simulation parameters

    # size of particle cube
    nParticlesX = 20  # number of particles in x direction
    nParticlesY = 20  # number of particles in y direction

    # distance between Particles at start
    distance = 0.5  # m

    # kernel radius
    h = np.sqrt(2) * distance  # m
    # number of particles
    N = nParticlesX * nParticlesY

    densityReference = 700  # kg/m^2

    nu = 0.2

    # mass particle
    mass = densityReference * distance**2  # kg

    # initial positions of particles
    pos = getPositionInGrid(N, nParticlesX, nParticlesY, distance)  # m
    # initial acceleration of particles
    acc = np.zeros((N, 2))  # m/s**2
    # initial velocity of particles
    vel = np.zeros((N, 2))  # m/s

    # show window
    fig, ax = plt.subplots()

    maxT = 200
    dt = 0.06

    oldT = 0
    oldTime = time.time()
    allDTs = []

    plt.show(block=False)

    # main loop
    for t in range(maxT):
        # (1/2) kick
        vel += acc * dt / 2

        # drift
        pos += vel * dt

        # collisions
        boundaryLeft = -8
        boundaryRight = 8
        boundaryTop = 1000
        boundaryBottom = 0

        boundDamping = -0.4

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
            s=40,
            alpha=0.5,
        )

        ax.set(xlim=(-60, 60), ylim=(-10, 110))
        plt.pause(0.001)

        if oldT < np.round((t / maxT) * 100):
            oldT = np.round((t / maxT) * 100)
            realDT = time.time() - oldTime
            oldTime = time.time()
            allDTs.append(realDT)

            print(
                oldT,
                "%, remaining:",
                np.round((realDT * (100 - oldT)) / 60, 3),
                "min; average dt:",
                sum(allDTs) / len(allDTs),
            )

    return 0


if __name__ == "__main__":
    main()
