import numpy as np
import matplotlib.pyplot as plt


N = 500
h = 100
m = 1

pos = np.random.randn(N,3)*100
p = np.zeros(N)

def W(r, h):
    return (1 / ((h**3) * (np.pi**(2/3)))) * np.exp((-(np.abs(r)**2)/h**2))

r_values = np.linspace(0, 100, 100)
W_values = W(r_values,h)

plt.plot(r_values, W_values)
plt.show()
