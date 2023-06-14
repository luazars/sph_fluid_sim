import numpy as np
import matplotlib.pyplot as plt


# setting positions of particles to a random position on the screen
def setRandomPosition(N):
    return np.random.random((N, 2))


# set position of particles in a cube
def getPositionInGrid(N, nParticlesX, nParticlesY, distance):
    # top-left position of particle cube
    xPos = 1 / 2 * (nParticlesX * distance)  # center the cube on x-axis
    yPos = 0.1

    pos = np.zeros((N, 2))
    pos = np.zeros((N, 2))
    x = np.tile(np.arange(nParticlesX), nParticlesY)
    y = np.repeat(np.arange(nParticlesY), nParticlesX)

    pos[:, 0] = x * distance + xPos
    pos[:, 1] = y * distance + yPos

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
    acc[:, 1] += -9.81
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

    ax.clear()
    ax.scatter(
        pos[:, 0],
        pos[:, 1],
        s=40,
        c=density,
        picker=True,
    )

    fig.canvas.draw()


def onKeyPress(
    event, ax, fig, acc, vel, pos, N, mass, density, densityReference, pressure, h
):
    if event.key == " ":
        timeStep(
            N, mass, 1, acc, vel, pos, density, densityReference, pressure, h, ax, fig
        )


def onPick(event, ax, fig, pos, density, pressure, acc, vel):
    ind = event.ind[0]  # Get the index of the clicked particle

    ax.clear()
    ax.scatter(
        pos[:, 0],
        pos[:, 1],
        s=40,
        c=density,
        picker=True,
    )

    x = pos[ind, 0]
    y = pos[ind, 1]

    # acc arrow
    a_x = acc[ind, 0]
    a_y = acc[ind, 1]
    ax.arrow(x, y, a_x, a_y, head_width=0.05, head_length=0.1, fc="red", ec="red")

    # vel arrow
    v_x = vel[ind, 0]
    v_y = vel[ind, 1]
    ax.arrow(x, y, v_x, v_y, head_width=0.05, head_length=0.1, fc="blue", ec="blue")

    fig.canvas.draw()

    ind = event.ind
    print("index:", ind, "pressure:", pressure[ind], "density:", density[ind])


def main():
    # Simulation parameters

    # size of particle cube
    nParticlesX = 20  # number of particles in x direction
    nParticlesY = 10  # number of particles in y direction

    # distance between Particles at start
    distance = 0.4  # m

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

    x = pos[:, 0]
    y = pos[:, 1]

    ax.scatter(
        x,
        y,
        s=40,
        c=density,
        picker=True,
    )

    for i in range(N):
        x = pos[i, 0]
        y = pos[i, 1]

        # acc arrow
        a_x = acc[i, 0]
        a_y = acc[i, 1]
        ax.arrow(x, y, a_x, a_y, head_width=0.05, head_length=0.1, fc="red", ec="red")

    fig.canvas.mpl_connect(
        "pick_event",
        lambda event: onPick(event, ax, fig, pos, density, pressure, acc, vel),
    )
    fig.canvas.mpl_connect(
        "key_press_event",
        lambda event: onKeyPress(
            event,
            ax,
            fig,
            acc,
            vel,
            pos,
            N,
            mass,
            density,
            densityReference,
            pressure,
            h,
        ),
    )

    plt.show()


if __name__ == "__main__":
    main()
