import numpy as np
import matplotlib.pyplot as plt
import random

# size of particle cube
nParticlesX = 20  # number of particles in x direction
nParticlesY = 20  # number of particles in y direction

# distance between Particles at start
distance = 0.4  # m

# kernel radiuss
h = np.sqrt(2) * distance  # m

# number of particles
N = nParticlesX * nParticlesY


# speed of sound
c = 5  # m/s
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
# velocity of particles
vel = np.zeros((N, 2))  # m/s
# acceleration of particles
acc = np.zeros((N, 2))  # m/s**2


# setting positions of particles to a random position on the screen
def setRandomPosition():
    for i in range(N):
        pos[i][0] = random.random()
        pos[i][1] = random.random()


# set position of particles in a cube
def setPositionInGrid():
    # top-left position of particle cube
    xPos = 1 / 2 * (nParticlesX * distance)  # center the cube on x-axis
    yPos = 0.1

    for i in range(N):
        x = i % nParticlesX
        y = int(i / N * nParticlesY)

        pos[i][0] = x * distance + xPos
        pos[i][1] = y * distance + yPos


# kernel function / W function
def W(r):
    # kernel radius, include all surrounding particles of one particle in a grid
    return (1.0 / (h * np.sqrt(np.pi))) ** 3 * np.exp(-(r**2) / h**2)


# gradient of the kernel function
def gradW(x, y, r):
    n = -2 * np.exp(-(r**2) / h**2) / h**5 / (np.pi) ** (3 / 2)
    wx = n * x
    wy = n * y
    return wx, wy


# calculate density of all particles
def calculateDensity():
    for i in range(N):
        for j in range(N):
            r = np.linalg.norm(pos[i] - pos[j])
            density[i] += mass * W(r)


# calculate pressure with density of all particles with tait equation
def calculatePressure():
    for i in range(N):
        pressure[i] = (
            (densityReference * c**2) / m * ((density[i] / densityReference) ** m - 1)
        )  # Pa


# calculate acceleration of all particles
def calculateAcceleration():
    calculateDensity()
    calculatePressure()

    for i in range(N):
        acc[i] = np.zeros(2)  # initialize acceleration for particle i

        for j in range(N):
            if i != j:
                r = pos[i] - pos[j]
                dist = np.linalg.norm(r)

                # Calculate the gradient of the kernel function
                grad = np.array(gradW(r[0], r[1], dist))

                # Calculate the pressure term contribution to acceleration
                pressure_term = (
                    -mass
                    * (
                        pressure[i] / (density[i] ** 2)
                        + pressure[j] / (density[j] ** 2)
                    )
                    * grad
                )

                # Accumulate acceleration due to pressure
                acc[i] += pressure_term


def showWindow():
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
        ax.arrow(
            pos[i, 0],
            pos[i, 1],
            acc[i, 0],
            acc[i, 1],
            head_width=0.05,
            head_length=0.1,
            fc="red",
            ec="red",
        )
    fig.canvas.mpl_connect("pick_event", lambda event: onPick(event, ax, fig))
    fig.canvas.mpl_connect("key_press_event", lambda event: onKeyPress(event, ax, fig))

    plt.show()


def onKeyPress(event, ax, fig):
    global vel, pos, acc
    dt = 1
    if event.key == " ":
        # (1/2) kick
        vel += acc * dt / 2

        # drift
        pos += vel * dt

        # update accelerations
        calculateAcceleration()

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


def onPick(event, ax, fig):
    ind = event.ind[0]  # Get the index of the clicked particle
    x = pos[ind, 0]
    y = pos[ind, 1]
    a_x = acc[ind, 0]
    a_y = acc[ind, 1]

    ax.clear()
    ax.scatter(
        pos[:, 0],
        pos[:, 1],
        s=40,
        c=density,
        picker=True,
    )

    ax.arrow(x, y, a_x, a_y, head_width=0.05, head_length=0.1, fc="red", ec="red")
    fig.canvas.draw()

    ind = event.ind
    print("index:", ind, "pressure:", pressure[ind], "density:", density[ind])


def main():
    setPositionInGrid()
    showWindow()


if __name__ == "__main__":
    main()
