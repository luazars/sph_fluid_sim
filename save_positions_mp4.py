import matplotlib.animation as animation
import time
import numpy as np
import matplotlib.pyplot as plt

positions = np.load("saved_positions/positions_20x200_nu.6.npy")
num_frames = positions.shape[0]

fig = plt.figure(figsize=(5, 5), dpi=(1920 / 5), facecolor="black")
ax = plt.gca()
ax.set_facecolor("#131621")


def update_pos(frame):
    ax.clear()
    ax.set(xlim=(-30.5, 30.5), ylim=(-0.5, 60))
    ax.set_aspect("equal")
    return ax.scatter(
        positions[frame, :, 0],
        positions[frame, :, 1],
        s=1,
        alpha=1,
        color="#99faff",
    )


simulation = animation.FuncAnimation(fig, update_pos, frames=num_frames, interval=100)
plt.close()

Writer = animation.writers["ffmpeg"]
writer = Writer(fps=30, metadata=dict(artist="Me"), bitrate=1800)
simulation.save("sph.mp4", writer=writer)
