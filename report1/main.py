import numpy as np
from numpy import typing as npt
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt
from enum import Enum, auto

type Vector = npt.ArrayLike[npt.float64]

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6


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
    

def Initalize(state: State) -> List[Atom]:
    """
        Returns a list of two atoms seperated by x Ångström.
        Depending on state the atoms are seperated by 4 or 3 Ångström.
    """
    if state == State.far:
        atoms: List[Atom] = [
            Atom(
                np.array([0,0,0]),
                np.array([0,0,0])
            ),
            Atom(
                np.array([4,0,0]),
                np.array([0,0,0])
         )
        ]
    else:
        atoms: List[Atom] = [
            Atom(
                np.array([0,0,0]),
                np.array([0,0,0])
            ),
            Atom(
                np.array([3,0,0]),
                np.array([0,0,0])
         )
        ]

    return atoms


def VerletMethodPosition(atom: Atom, dt: float, force: Vector) -> None:
    """
        Computes the new position of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit eV/c^2
    """
    pos: Vector = atom.position + atom.velocity * dt + force / (2 * atom.mass) * dt**2

    atom.position = pos


def VerletMethodVelocity(atom: Atom, dt: float, forceNew: Vector, forceOld: Vector) -> None:
    """
        Computes the new velocity of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit eV/c^2
    """
    vel: Vector = atom.velocity + (forceNew  + forceOld) / (2 * atom.mass) * dt
    atom.velocity = vel


def KineticEnergy(atom: Atom, velocity: Vector) -> float:
    """
        Computes the kinetic energy of an atom.
    """
    return 0.5 * atom.mass * np.linalg.norm(velocity)**2


def computeHamiltonian(atoms: List[Atom], velocity: Vector, potential: float) -> float:
    """
        Computes the Hamiltonian of the system.
    """
    kinetic_part: float = 0
    for i in range(len(atoms)):
        kinetic_part += atoms[i].mass * np.linalg.norm(velocity[i])**2
    
    return 1/2 * kinetic_part + potential


def plotSimulation(t: np.ndarray, atom1_pos_x: List[float], atom2_pos_x: List[float],
                atom1_vel_x: List[float], atom2_vel_x: List[float], atom1_Kinetic: List[float],
                atom2_Kinetic: List[float], Potential: List[float], Hamiltonian: List[float]) -> None:
    """
        Function that plots the different quantities in the system.
        Then stores the plots as .png files.
    """

    plt.title('Position of the atoms')
    plt.plot(t, atom1_pos_x, label='$x_1(t)$')
    plt.plot(t, atom2_pos_x, label='$x_2(t)$')
    plt.xlabel('Time [s]')
    plt.ylabel('Position [Å]')
    plt.legend()
    plt.grid()
    plt.savefig('Position of the atoms.png')
    #plt.show()

    plt.clf()

    plt.title('Velocity of the atoms')
    plt.plot(t, atom1_vel_x, label='$v_1(t$)')
    plt.plot(t, atom2_vel_x, label='$v_2(t)$')
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [Å/s]')
    plt.legend()
    plt.grid()
    plt.savefig('Velocity of the atoms.png')
    # plt.show()

    plt.clf()

    plt.title('Different energies in the system')

    plt.plot(t, np.array(atom2_Kinetic) + np.array(atom1_Kinetic), label='$\mathcal{K}(t)$')

    plt.plot(t, Hamiltonian, label='$\mathcal{H}(t)$')

    plt.plot(t, Potential, '--',  markersize = 1, label='$\mathcal{V}(t)$')

    plt.grid()
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [eV]')
    plt.legend()
    plt.savefig('Different energies in the system.png')
    # plt.show()


def runSimulation(plotSim: bool, plotHamiltonian: bool, state: State, dt = 1e-15):
    """
        Run the simulation.
        plotSim: determines if one wants to plot the energy, position and velocity plot
        plotHamiltonian: recurise call to plot the hamiltonian for various dts
        state: used to determine which position to seperate the atoms about.
    """
    atoms = Initalize(state)

    T: float = 1/1e11
    t: float = 0

    if plotHamiltonian:
        """

        """
        dt = [1e-12, 5e-13, 1e-13, 5e-14]
        for i in range(len(dt)):
            Hamiltonian = runSimulation(False, False, state, dt[i])
            plt.plot(np.linspace(0, T, len(Hamiltonian)), Hamiltonian, label=f'dt = ${dt[i]}$')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [eV]')
        plt.title('Hamiltonian with different time-steps')
        plt.grid()

        plt.savefig('Hamiltonian with different time-steps.png')

        return
        
    
    # Various lists that stores the elements from the simulation
    Potential: List[float] = []

    atom1_pos_x: List[float] = []
    atom2_pos_x: List[float] = []

    atom1_vel_x: List[float] = []
    atom2_vel_x: List[float] = []

    atom1_Kinetic: List[float] = []

    atom2_Kinetic: List[float] = []

    Hamiltonian: List[float] = []

    """
        While loop that runs the simulation. 
        Iterates via the Velocity Verlet method for two atoms.
    """
    atom1: Atom = atoms[0]
    atom2: Atom = atoms[1]
    while t < T:
        """
            Loop that iterates over the Velocity Verlet method
        """

        potential, force1 = LennardJones(atom1.position, atom2.position)
        force2 = - force1 # Newtons third law
        print(force1)

        VerletMethodPosition(atom1, dt, force1)
        VerletMethodPosition(atom2, dt, force2)
        
        atom1_pos_x.append(atom1.position[0])
        atom2_pos_x.append(atom2.position[0])

        _, newForce1 = LennardJones(atom1.position, atom2.position)
        newForce2 = - newForce1 # Newtons third law

        atom1_vel_x.append((atom1.velocity[0]))
        atom2_vel_x.append((atom2.velocity[0]))

        VerletMethodVelocity(atom1, dt, newForce1, force1)
        VerletMethodVelocity(atom2, dt, newForce2, force2)

        print(atom1.velocity)

        Potential.append(potential)

        atom1_Kinetic.append(KineticEnergy(atom1, atom1.velocity))
        atom2_Kinetic.append(KineticEnergy(atom2, atom2.velocity))

        Hamiltonian.append(computeHamiltonian(atoms, [atom1.velocity, atom2.velocity], potential))


        print(f'Completed {t/T} %')
            
        t += dt

    t: Vector = np.linspace(0, T, len(atom1_pos_x)) # Time array

    if plotSim:
        plotSimulation(t, atom1_pos_x, atom2_pos_x, atom1_vel_x, atom2_vel_x, atom1_Kinetic, atom2_Kinetic, Potential, Hamiltonian)

    return Hamiltonian


if __name__ == '__main__':
    computeForceAndPotentialForPlot()
    state = State.far
    runSimulation(False, True, state)
    runSimulation(True, False, state)
    state = State.close
    # runSimulation(True, False, state) # Used for investigating atoms placed close together
