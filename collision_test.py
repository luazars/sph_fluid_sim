import numpy as np
import matplotlib.pyplot as plt


def getPositionInGrid(N, nParticlesX, nParticlesY, distance):
    # top-left position of particle cube
    xPos = (nParticlesX * distance) / 2  # center the cube on x-axis

    pos = np.zeros((N, 2))

    x = np.tile(np.arange(nParticlesX), nParticlesY)
    y = np.repeat(np.arange(nParticlesY), nParticlesX)

    pos[:, 0] = x * distance - xPos
    pos[:, 1] = y * distance

    return pos


def main():
    py = 10
    px = 10
    N = py * px
    # initial positions of particles
    pos = getPositionInGrid(N, py, px, 2)  # m
    # initial acceleration of particles
    acc = np.zeros((N, 2))  # m/s**2
    # initial velocity of particles
    vel = np.random.random((N, 2))  # m/s

    # show window
    fig, ax = plt.subplots()

    maxT = 1000
    dt = 0.1

    # main loop
    for t in range(maxT):
        acc[:][1] = -1
        vel += acc * dt
        pos += vel * dt

        # collisions
        boundaryLeft = -10
        boundaryRight = 10
        boundaryTop = 100
        boundaryBottom = -3

        if t > 400:
            boundaryRight = 6
            boundaryLeft = -6
        if t > 800:
            boundaryBottom = -8

        boundDamping = -0.5
        boundFix = 0.1

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

        # update window
        plt.sca(ax)
        plt.cla()
        plt.scatter(
            pos[:, 0],
            pos[:, 1],
            s=40,
            alpha=0.5,
        )
        ax.set(xlim=(-20, 20), ylim=(-10, 20))
        plt.pause(0.0001)

    return 0


if __name__ == "__main__":
    main()
