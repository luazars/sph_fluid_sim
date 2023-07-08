import numpy as np
import matplotlib.pyplot as plt


def main():
    positions = np.load("saved_positions/positions_20x200_nu.6.npy")

    num_frames = positions.shape[0]
    _, ax = plt.subplots(facecolor="#000000")
    ax.set_facecolor("#131621")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_color("black")

    for frame in range(num_frames):
        ax.clear()
        ax.set(xlim=(-30, 30), ylim=(-1, 60))
        plt.gca().set_aspect("equal")
        ax.scatter(
            positions[frame, :, 0],
            positions[frame, :, 1],
            s=1,
            alpha=1,
            color="#99faff",
        )
        plt.pause(0.01)

    plt.show()
    return 0


if __name__ == "__main__":
    main()
