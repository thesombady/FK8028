import numpy as np
from numpy import typing as npt
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt
from enum import Enum, auto

type Vector = npt.ArrayLike[npt.float64]

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6

IMAGESIZE = 10.0 # Ångström
DENSITY = 0.001 # Ångström^-3

class State(Enum):
    """
        Enum to choose starting relative position
    """
    close = auto()
    far = auto()
    

@dataclass
class Atom: # Argon atom
    position: Vector
    velocity: Vector 
    radius: float = 1.0 # Ångström
    mass: float = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2


def LennardJones(pos1: Vector, pos2: Vector) -> Tuple[float, np.ndarray]:
    """
       Computes the Lennards-Jones potential and force between two atoms. 
    """
    sigma: float = 3.40 # Å
    epsilon: float = 0.0104 # eV

    r: Vector = pos1 - pos2 # Relative position

    delta: float = (sigma / np.linalg.norm(r)) ** 6
    
    potential: float = 4 * epsilon * (delta ** 2 - delta)

    force: Vector = 4 * epsilon * (12 * delta ** 2 - 6 * delta) * (r) / np.linalg.norm(r)**2

    return potential, force

    
def plotForcePotential(V: List[float], F_x: List[float], range_: Vector) -> None:
    """
        Plots the force and potential given a potential and force in the X-direction.
        Also plots the force as a differece quotient.
    """
    plt.title('Interatomic force in the x direction')
    plt.plot(range_, F_x, label='${F_x}$')
    plt.xlabel('Distance [Å]')
    plt.ylabel('Force [eV/Å]')
    plt.plot(range_[:-1], [(V[i+1]-V[i])/(range_[1]-range_[0]) for i in range(len(V) - 1)], 
            '.', markersize = 2, label='${F_x^{\Delta}}$'
    )
    plt.grid()
    plt.legend()
    plt.savefig('Interatomic force in the x direction.png')
    plt.clf()
    # plt.show()

    plt.title('Lennard-Jones potential')
    plt.plot(range_, V)
    plt.grid()
    plt.xlabel('Distance [$Å$]')
    plt.ylabel('Potential [eV]')
    plt.savefig('Lennard-Jones potential.png')
    # plt.show()


def computeForceAndPotentialForPlot() -> None:
    """
        Computes the force in the x direction for each atom in the system.
        however only the x component is computed here.
    """ 
    Atom1 = np.array([0,0,0])
    
    fx: List[float] = []
    V: List[float] = []
    range_: np.ndarray = np.linspace(3, 10, 150)

    for i in range_:
        Atom2 = np.array([i, 0, 0])
        potential, force = LennardJones(Atom1, Atom2) # Computes the lennard jones potential and force
        fx.append(force[0])
        V.append(potential)
    
    plotForcePotential(V, fx, range_)
    

def Initalize(n: int) -> List[Atom]:
    """
        Returns a list of two atoms seperated by x Ångström.
        Depending on state the atoms are seperated by 4 or 3 Ångström.
    """


    return [
        Atom(
            np.array([0,0,0]),
            np.array([0,0,0])
        ),
        Atom(
            np.array([4,0,0]),
            np.array([0,0,0])
        )
    ]


def VerletMethodPosition(atom: Atom, dt: float, acc: Vector) -> Vector:
    """
        Computes the new position of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit eV/c^2
    """
    pos: Vector = atom.position + atom.velocity * dt + acc / (2) * dt**2

    return pos


def VerletMethodVelocity(atom: Atom, dt: float, accNew: Vector, accOld: Vector) -> Vector:
    """
        Computes the new velocity of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit eV/c^2
    """
    vel: Vector = atom.velocity + (accNew  + accOld) / (2) * dt
    return vel


def validateBoundary(pos: Vector) -> Vector:
    """
        Function that checks if the position of an atom is within the boundary.
        If not it will return the position of the atom on the other side of the boundary.
    """

    pos_copy = pos.copy()
    for i in range(3):
        if pos[i] > IMAGESIZE:
            pos_copy[i] = pos[i] - IMAGESIZE
        elif pos[i] < 0:
            pos_copy[i] = IMAGESIZE + pos[i]
    
    return pos_copy 


def runSimulation(numberOfAtoms: int = 2):
    """
        Run the simulation.
        plotSim: determines if one wants to plot the energy, position and velocity plot
        plotHamiltonian: recurise call to plot the hamiltonian for various dts
        state: used to determine which position to seperate the atoms about.
    """
    atoms = Initalize(numberOfAtoms)
    T = 1e-11
    t = 0
    dt = 1e-15

    """
        While loop that runs the simulation. 
        Iterates via the Velocity Verlet method for two atoms.
    """

    Accelerations = []

    # Computes the initial acceleration for each atom to avoid recomputing it in the loop

    for i in range(len(atoms)):
        force: Vector = np.array([0.0, 0.0, 0.0])
        for j in range(len(atoms)):
            if i == j:
                continue
            _, force_ = LennardJones(atoms[i].position, atoms[j].position)
            force += force_
        
        Accelerations.append(force / atoms[i].mass)
    
    pos1_x: List[float] = []
    pos2_x: List[float] = []

    while t < T:
        """
            Loop that iterates over the Velocity Verlet method
        """
        
        pos1_x.append(atoms[0].position[0])
        pos2_x.append(atoms[1].position[0])

        
        for i in range(len(atoms)):
            pos: Vector = VerletMethodPosition(atoms[i], dt, Accelerations[i])
            pos = validateBoundary(pos)
            atoms[i].position = pos 
            

        
        newAccelerations = []

        # Do a map that contains the (i,j) combination of atoms that have already been computed

        for i in range(len(atoms)):
            newForce: Vector = np.array([0.0, 0.0, 0.0])
            for j in range(len(atoms)):
                if i == j:
                    continue
                _, force_ = LennardJones(atoms[i].position, atoms[j].position)
                newForce += force_
            
            print(newForce) # [0.00484532 0.         0.        ]

            atoms[i].velocity = VerletMethodVelocity(atoms[i], dt, newForce/atoms[i].mass, Accelerations[i])

            print(atoms[i].velocity) # [3.31772493e+11 0.00000000e+00 0.00000000e+00] # Same order of magnitude as previous simulation

            newAccelerations.append(newForce / atoms[i].mass)

        Accelerations = newAccelerations

        t += dt

    plt.plot(np.linspace(0, T, len(pos1_x)), pos1_x, label='atom 1')
    plt.plot(np.linspace(0, T, len(pos2_x)), pos2_x, label = 'atom2')
    plt.legend()
    plt.show()
    


if __name__ == '__main__':
    state = State.far
    runSimulation()
