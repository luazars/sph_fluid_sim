import numpy as np

# Create an array of densities
densities = np.array([1000, 1100, 1200])

def tait_equation(density, B=0.3214, n=7.0, P0=0.0, rho0=1.0):
    """
    Calculates the pressure of a fluid based on its density using the Tait equation.

    Parameters:
    density (float or numpy array): the density of the fluid
    B (float): a constant for the fluid (default value for water: 0.3214)
    n (float): a constant for the fluid (default value for water: 7.0)
    P0 (float): the reference pressure (default: 0.0)
    rho0 (float): the reference density (default: 1.0)

    Returns:
    The pressure of the fluid.
    """
    return P0 + B * ((density / rho0) ** n - 1)

# Calculate the pressure using the Tait equation
pressures = tait_equation(densities)

# Print the results
print(pressures)
