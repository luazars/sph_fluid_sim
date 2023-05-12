import numpy as np
import matplotlib.pyplot as plt

h = 0.05 * np.sqrt(2)


def W(r, h):
    return (1 / ((h**3) * (np.pi ** (2 / 3)))) * np.exp((-(np.abs(r) ** 2) / h**2))


r_values = np.linspace(0, 0.1, 100)
W_values = W(r_values, h)

plt.plot(r_values, W_values / W(0, h))

plt.show()
