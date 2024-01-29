import numpy as np
import matplotlib.pyplot as plt


# Constants for Argon
sigma = 3.4e-10  # m
epsilon = 1.65e-21  # J
mass = 6.69e-26  # kg
timestep = 1e-14  # s

def lennard_jones(r1, r2):
    """
    Calculate the Lennard-Jones potential and force for a given distance.
    """
    r = r1 - r2
    fraction = (sigma/r) ** 6
    V = 4*epsilon * ((fraction)**2 - (fraction))
    F = 4*epsilon* r/ r**2 * (12 * fraction**2 - 6* fraction)
    return V, F

# Initial conditions
r1 = 0
v1 = 0

r2 = 4e-10
v2 = 0
_, a = lennard_jones(r1, r2)

a = a/mass

def velocity_verlet_pos(r, v, a):
    """
    Perform one step of the Velocity Verlet algorithm with the first atom fixed.
    """
    return r + v*timestep + 0.5*a*timestep**2

def velocity_verlet_vel(v, aNew, aOld):
    """
    Perform one step of the Velocity Verlet algorithm with the first atom fixed.
    """
    return v + 0.5*(aNew + aOld)*timestep


# Perform the Velocity Verlet algorithm for 100 steps

T = 300 * timestep
t = 0
x1 = []
x2 = []

p1 = []
p2 = []

H = []
P = []
while t < T:
    p, a = lennard_jones(r1, r2) 
    a = a/mass
    r1_new = velocity_verlet_pos(r1, v1, a)
    r2_new = velocity_verlet_pos(r2, v2, -a)

    _, aNew = lennard_jones(r1_new, r2_new)
    aNew = aNew/mass

    H.append(mass*(v1**2 + v2**2) + p)
    P.append(p)

    v1 = velocity_verlet_vel(v1, aNew, a)
    v2 = velocity_verlet_vel(v2, -aNew, -a)
    print(P)

    x1.append(r1)
    x2.append(r2)

    p1.append(v1)
    p2.append(v2)

    r1 = r1_new
    r2 = r2_new

    a = aNew

    t += timestep
plt.plot(np.linspace(0, T, len(x1)), x1, label='Atom 1')
plt.plot(np.linspace(0, T, len(x2)), x2, label='Atom 2')
plt.show()

plt.plot(np.linspace(0, T, len(p1)), H, label='Hamiltonian')
plt.plot(np.linspace(0, T, len(p1)), P, label='Potential')
plt.legend()
plt.show()
