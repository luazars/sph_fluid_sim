import numpy as np
import matplotlib.pyplot as plt


def main():
    positions = np.load("saved_positions/positions_20x200_nu.6.npy")

    num_frames = positions.shape[0]
    _, ax = plt.subplots()

    for frame in range(num_frames):
        ax.clear()
        ax.set(xlim=(-40, 40), ylim=(0, 80))
        plt.gca().set_aspect("equal")
        ax.scatter(positions[frame, :, 0], positions[frame, :, 1], s=5, alpha=0.5)
        plt.pause(0.01)

    plt.show()
    return 0


if __name__ == "__main__":
    main()
