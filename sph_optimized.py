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
    pos[:, 1] = y * distance + 10

    colors = np.full(N, "black")

    # for i in range(N):
    #     if pos[i, 0] < 0:
    #         pos[i, 0] -= 3
    #         if pos[i, 1] > 10:
    #             pos[i, 1] += 10
    #             colors[i] = "red"
    #         else:
    #             colors[i] = "blue"

    #     else:
    #         pos[i, 0] += 3
    #         if pos[i, 1] > 10:
    #             pos[i, 1] += 10
    #             colors[i] = "yellow"
    #         else:
    #             colors[i] = "green"

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
def getPressure(density, densityReference, pressureReference):
    c = 5  # speed of sound
    m = 7  # Adiabatenexponent
    B = (densityReference * c**2) / m

    pressure = B * (((density / densityReference) ** m) - 1) + pressureReference
    pressure = np.maximum(pressure, 0)

    return pressure  # Pa


# calculate acceleration of all particles
# https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1
def getAcceleration(N, mass, pos, vel, densityReference, pressureReference, h, nu, g):
    density = getDensity(mass, N, pos, h)
    pressure = getPressure(density, densityReference, pressureReference)

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
    acc += g

    # viscosity
    acc -= nu * vel

    return acc


def main():
    # Simulation parameters
    simulateInRealtime = False

    # size of particle cube
    nParticlesX = 20  # number of particles in x direction
    nParticlesY = 200  # number of particles in y direction

    # distance between Particles at start
    distance = 0.5  # m

    densityReference = 700  # kg/m^2
    pressureReference = 0

    #  weight force
    g = (0, -2)

    # viscosity
    nu = 0.2

    maxT = 400
    dt = 0.06

    # number of particles
    N = nParticlesX * nParticlesY
    # kernel radius
    h = np.sqrt(2) * distance  # m
    # mass particle
    mass = densityReference * distance**2  # kg

    # initial positions of particles
    pos, colors = getPositionAndColorInGrid(N, nParticlesX, nParticlesY, distance)  # m
    # initial velocity of particles
    vel = np.zeros((N, 2))  # m/s
    # initial acceleration of particles
    acc = np.zeros((N, 2))  # m/s**2

    # show window
    fig, ax = plt.subplots()

    if simulateInRealtime == False:
        positions = []
        oldT = 0
        oldTime = time.time()
        allDTs = []
    else:
        plt.show(block=False)

    # main loop
    for t in range(maxT):
        # (1/2) kick
        vel += acc * dt / 2

        # drift
        pos += vel * dt

        # collisions
        boundaryLeft = -40
        boundaryRight = 40
        boundaryTop = 1000
        boundaryBottom = 0

        boundDamping = -0.5
        boundFix = 0.5

        # Check if particles collide with the boundaries
        out_of_bounds_left = pos[:, 0] <= boundaryLeft
        out_of_bounds_right = pos[:, 0] >= boundaryRight
        out_of_bounds_top = pos[:, 1] >= boundaryTop
        out_of_bounds_bottom = pos[:, 1] <= boundaryBottom

        # Reflect velocities
        vel[out_of_bounds_left, 0] *= boundDamping
        vel[out_of_bounds_right, 0] *= boundDamping
        vel[out_of_bounds_top, 1] *= boundDamping
        vel[out_of_bounds_bottom, 1] *= boundDamping

        # Add small vel to prevent particles of getting stuck in the border
        vel[out_of_bounds_left, 0] += boundFix
        vel[out_of_bounds_right, 0] -= boundFix
        vel[out_of_bounds_top, 1] -= boundFix
        vel[out_of_bounds_bottom, 1] += boundFix

        # Update positions
        pos[out_of_bounds_left, 0] = boundaryLeft
        pos[out_of_bounds_right, 0] = boundaryRight
        pos[out_of_bounds_top, 1] = boundaryTop
        pos[out_of_bounds_bottom, 1] = boundaryBottom

        # update accelerations
        acc = getAcceleration(
            N, mass, pos, vel, densityReference, pressureReference, h, nu, g
        )

        # (1/2) kick
        vel += acc * dt / 2

        if simulateInRealtime:
            # update window
            plt.sca(ax)
            plt.cla()
            plt.scatter(
                pos[:, 0],
                pos[:, 1],
                s=40,
                alpha=0.5,
                c=colors,
            )
            ax.set(xlim=(-40, 40), ylim=(0, 80))
            plt.gca().set_aspect("equal")
            plt.pause(0.001)
        else:
            positions.append(pos.copy())

            # print rendering progress
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

    if simulateInRealtime == False:
        # save positions
        positions = np.array(positions)
        np.save("saved_positions/positions.npy", positions)

    return 0


if __name__ == "__main__":
    main()
