import numpy as np
import matplotlib.pyplot as plt



N = 500
h = 1
m = 1

pos = np.random.randn(N,3)*100
p = np.zeros(N)

def W(r, h):
    return (1 / (h * np.sqrt(2*np.pi))) * np.exp((-(r**2)/2*h**2))

def R(x,y,z):
    return np.sqrt(x**2+y**2+z**2)

def V(h):
    return (4.0/3.0)*np.pi*h**3


# Define a range of values for r
r_values = np.linspace(0, 100, 100)

# Choose a value for h
h = 0.1

# Calculate the W function for each value of r
W_values = W(r_values,h)

# Plot the results
plt.plot(r_values, W_values)
plt.show()
