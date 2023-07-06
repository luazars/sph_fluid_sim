import numpy as np
import matplotlib.pyplot as plt


def main():
    positions = np.load("saved_positions/positions_50x100_nu.1.npy")

    num_frames = positions.shape[0]
    _, ax = plt.subplots()

    for frame in range(num_frames):
        ax.clear()
        ax.set(xlim=(-50, 50), ylim=(-10, 50))
        ax.scatter(positions[frame, :, 0], positions[frame, :, 1], s=30, alpha=0.5)
        plt.pause(0.000000000001)
        print(frame)

    plt.show()
    return 0


if __name__ == "__main__":
    main()
