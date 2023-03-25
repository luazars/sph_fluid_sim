import numpy as np
import matplotlib.pyplot as plt


N = 500
h = .1
m = 1

pos = np.random.randn(N,3)*100
p = np.zeros(N)

def W(r, h):
    return (1 / (h * np.sqrt(2*np.pi))) * np.exp((-(r**2)/2*h**2))

r_values = np.linspace(0, 100, 100)
W_values = W(r_values,h)

plt.plot(r_values, W_values)
plt.show()
