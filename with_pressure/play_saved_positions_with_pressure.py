import numpy as np
import matplotlib.pyplot as plt


def main():
    positions = np.load("saved_positions/positions.npy")
    pressures = np.load("saved_pressures/pressures.npy")

    num_frames = positions.shape[0]

    _, ax = plt.subplots(facecolor="black")
    ax.set_facecolor("#131621")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")

    for frame in range(num_frames):
        ax.clear()
        ax.set(xlim=(-20.5, 40.5), ylim=(-0.5, 60.5))
        plt.gca().set_aspect("equal")
        ax.scatter(
            positions[frame, :, 0],
            positions[frame, :, 1],
            s=1,
            alpha=1,
            c=pressures[frame],
        )
        plt.pause(0.01)

    plt.show()
    return 0


if __name__ == "__main__":
    main()
