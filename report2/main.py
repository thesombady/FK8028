import numpy as np
from numpy import typing as npt
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt
from enum import Enum, auto

type Vector = npt.ArrayLike[npt.float64] # Vector type used when type-hinting

"""
    Constants used in the simulation.
"""

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6

GRAM_TO_EV = 5.6095886 * 1e32 / (SPEED_OF_LIGHT * 1e10)**2

DENSITY = 1.40 # g/cm^3

DENSITY = DENSITY * GRAM_TO_EV / (1e24) # eVs^2/Å^2 * 1/Å^3

BOXSIZE = 10.0 # Ångström

ATOM_MASS = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2

BOLTZMANN = 8.617333262145 * 1e-5 # eV/K


class State(Enum):
    """
       Enum for choosing which system to simulate.
       Bouncing test or lattice test.
    """
    bouncing = auto()
    lattice = auto()
    
    
@dataclass
class Atom: # Argon atom
    name: str
    position: Vector
    velocity: Vector 
    radius: float = 1.0 # Ångström
    mass: float = ATOM_MASS # eVs^2/Å^2


def LennardJones(pos1: Vector, pos2: Vector) -> Tuple[float, np.ndarray]:
    """
       Computes the Lennards-Jones potential and force between two atoms. 
    """
    sigma: float = 3.40 # Å
    epsilon: float = 0.0104 # eV

    r: Vector = pos1 - pos2 # Relative position

    for i in range(3):
        r[i] = minimumImageConvension1D(pos1[i], pos2[i])

    delta: float = (sigma / np.linalg.norm(r)) ** 6
    
    potential: float = 4 * epsilon * (delta ** 2 - delta)

    force: Vector = -4 * epsilon * (12 * delta ** 2 - 6 * delta) * (r) / np.linalg.norm(r)**2

    return potential, force


def Initalize(numberOfAtoms1D: int, state: State) -> List[Atom]:
    """
        : Returns a list of atoms in a cubic lattice.
        * numberOfAtoms: The number of atoms in one direction of the lattice in SC.
        
    """
    global BOXSIZE

    if state == State.bouncing:
        return [
            Atom(
                'Ar:0',
                np.array([0.07,0,0]),
                np.array([0,0,0])
            ),
            Atom(
                'Ar:1',
                np.array([6.07,0,0]),
                np.array([0,0,0])
            )
        ]

    BOXSIZE = ( ATOM_MASS * numberOfAtoms1D ** 3 / DENSITY ) **  ( 1/ 3)

    atoms: List[Atom] = []

    spacing = BOXSIZE / (numberOfAtoms1D - 1) * 0.85 # Spacing between the atoms

    print('Size: ' + str(BOXSIZE))

    counter = 0

    for i in range(numberOfAtoms1D):
        for j in range(numberOfAtoms1D):
            for k in range(numberOfAtoms1D):

                velocity = np.random.standard_normal(3) * 1e-10 # Å/s 'random' velocity

                atoms.append(
                    Atom(
                        f'Ar:{counter}',
                        np.array([1 + (i) * spacing, 1 + (j) * spacing, 1 + (k) * spacing]), # Position
                        velocity # Velocity
                    )
                )
                counter += 1
                
    return atoms



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


def minimumImageConvension1D(pos1: float, pos2: float) -> float:
    """
        Function that computes the minimum image convension for one dimension.
        * pos1: The position of the first atom in the any dimension
        * pos2: The position of the first atom in the any dimension
    """
    delta = pos2 - pos1
    delta += -BOXSIZE* round(delta / BOXSIZE)
    return delta


def validateBoundary(pos: Vector) -> Vector:
    """
        Function that checks if the position of an atom is within the boundary.
        If not it will return the position of the atom on the other side of the boundary.
    """
    
    return pos - np.floor(pos / BOXSIZE) * BOXSIZE 

def kineticTemperature(velocities: List[Vector]) -> Tuple[float, float]:
    """
        Computes the kinetic energy, and temperature of the system.
        * velocities: List of all the velocities of the atoms in the system.
    """
    K = 0
    
    for v in velocities:
        K += 0.5 * ATOM_MASS * v**2
    
    T = 2 * K / (3 * BOLTZMANN * len(velocities))

    return K, T


def runSimulation(state: State, numberOfAtoms: int = 3):
    """
        Run the simulation.
        plotSim: determines if one wants to plot the energy, position and velocity plot
        plotHamiltonian: recurise call to plot the hamiltonian for various dts
        state: used to determine which position to seperate the atoms about.
    """
    atoms = Initalize(numberOfAtoms, state)
    if state == State.bouncing:
        T = 1e-11
        dt = 1e-16
    else:
        T = 1e-12 # 1e-12
        dt = 1e-16
    t = 0

    """
        While loop that runs the simulation. 
        Iterates via the Velocity Verlet method for two atoms.
    """

    ax = plt.figure().add_subplot(projection='3d')
    ax.set_title('Initial position of atoms')

    for i in range(len(atoms)):
        ax.scatter(atoms[i].position[0], atoms[i].position[1], atoms[i].position[2])
    ax.set_xticks(np.arange(0, BOXSIZE, numberOfAtoms))
    ax.set_yticks(np.arange(0, BOXSIZE, numberOfAtoms))
    ax.set_zticks(np.arange(0, BOXSIZE, numberOfAtoms))
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_zlabel('z [Å]')
    plt.savefig('initial_position.png')
    plt.clf()

    Accelerations: Vector = []
    Potential: float = []

    # Computes the initial acceleration for each atom to avoid recomputing it in the loop

    Potential: List[float] = []

    allAtoms = {atom.name: {'position': [], 'velocity': [], 'potential': []} for atom in atoms}


    Potential_: List[float] = []
    for i in range(len(atoms)):
        
        name: str = 'Ar:{}'.format(i)

        force: Vector = np.array([0.0, 0.0, 0.0])
        
        potential: float = 0

        for j in range(len(atoms)):
            if i == j:
                continue
            p, force_ = LennardJones(atoms[i].position, atoms[j].position)
            force += force_
            potential += p

        (allAtoms[name]['potential']).append(potential)

        Accelerations.append(force / atoms[i].mass)
        Potential_.append(potential)

    Potential.append(sum(Potential_))

    Temperature: List[float] = []

    counter = 0
    while t < T:
        """
            Loop that iterates over the Velocity Verlet method
        """
        
        # Compute the new position for each atom
        for i in range(len(atoms)):
            name: str = 'Ar:{}'.format(i)
            pos: Vector = VerletMethodPosition(atoms[i], dt, Accelerations[i])
            pos = validateBoundary(pos)
            atoms[i].position = pos 
            (allAtoms[name]['position']).append(atoms[i].position)
        
        newAccelerations = []

        computedValues = dict()

        Potential_: List[float] = []
        # Compute the new acceleration for each atom, and thus also the new velocity
        for i in range(len(atoms)):
            name: str = 'Ar:{}'.format(i)
            
            potential: float = 0
            newForce: Vector = np.array([0.0, 0.0, 0.0])
            for j in range(len(atoms)):
                if i == j:
                    continue
                if computedValues.get((j, i)) is not None:
                    p, force_ = computedValues[(j, i)]
                    p = 0 # we don't want to double count the potential
                    force_ *= -1
                else:
                    p, force_ = LennardJones(atoms[i].position, atoms[j].position)
                    computedValues[(i, j)] = (p, force_)
                potential += p
                newForce += force_
            
            Potential_.append(potential)
            
            atoms[i].velocity = VerletMethodVelocity(atoms[i], dt, newForce/atoms[i].mass, Accelerations[i])

            (allAtoms[name]['velocity']).append(atoms[i].velocity)            

            newAccelerations.append(newForce / atoms[i].mass)

            (allAtoms[name]['potential']).append(potential)
        
        Potential.append(sum(Potential_))
        
        velo: List[Vector] = [allAtoms[name]['velocity'][-1] for name in allAtoms] # List of all current velocites

        K, Temp = kineticTemperature(velo)

        Temperature.append(Temp)
        
        # Rescale the velocities to the desired temperature
        if t < 1e-13 and counter % 50 == 0:

            for i in range(len(atoms)):
                velo[i] *= np.sqrt(94.4 / Temp)

            K, Temp = kineticTemperature(velo)
            Temperature[-1] = Temp

            for i in range(len(atoms)):
                name: str = 'Ar:{}'.format(i)
                atoms[i].velocity = velo[i]
                allAtoms[name]['velocity'][-1] = velo[i]
                

        Accelerations = newAccelerations

        t += dt
        counter += 1
        print(f'Completed {round(t/T, 4) * 100} %')

    t = np.linspace(0, T, len(allAtoms['Ar:0']['position']))

    if state == State.bouncing:
        """
            Plot the various plots for the Bouncing test.
        """
        # Plot of the position of the two atoms
        plt.title('Position of two particles')
        pos1 = np.array([allAtoms['Ar:0']['position'][i][0] for i in range(len(allAtoms['Ar:0']['position']))])
        pos2 = np.array([allAtoms['Ar:1']['position'][i][0] for i in range(len(allAtoms['Ar:1']['position']))])
        plt.plot(t, pos1, label="Atom 1")
        plt.plot(t, pos2, label="Atom 2")
        plt.xlabel("Time [s]")
        plt.ylabel("Position [Å]")
        plt.legend()
        plt.savefig('bsX.png')
        #plt.show()
        plt.clf()
        
        vel1 = np.array([allAtoms['Ar:0']['velocity'][i][0] for i in range(len(allAtoms['Ar:0']['position']))])
        vel2 = np.array([allAtoms['Ar:1']['velocity'][i][0] for i in range(len(allAtoms['Ar:1']['position']))])

        # Plot of the velocity of the two atoms
        plt.title('Velocity of two particles')
        plt.plot(t, vel1, label="Atom 1")
        plt.plot(t, vel2, label="Atom 1")
        plt.xlabel("Time [s]")
        plt.ylabel("Velocity [Å/s]")
        plt.legend()
        plt.savefig('bsV.png')
        #plt.show()
        plt.clf()

        Kinetic = atoms[0].mass * (vel1)**2+ atoms[1].mass * (vel2)**2

        potential = np.array([allAtoms['Ar:0']['potential'][j]for j in range(len(vel1))])
        # Plot of the energy of the two atoms
        plt.title('Energy of two particles')
        plt.plot(t, Kinetic, label="$\mathcal{K}(t)$")
        plt.plot(t, Kinetic + potential, label="$\mathcal{H}(t)$")
        plt.plot(t, potential, label="$\mathcal{V}(t)$")
        plt.xlabel("Time [s]")
        plt.ylabel("Velocity [eV]")
        plt.grid()
        plt.legend()
        plt.savefig('bsE.png')
        #plt.show()
        plt.clf()
        return 
        

    """
        Plot the various MD plots for the lattice test.
    """


    simulationLength = len(allAtoms['Ar:0']['position'])
    # Plot of the position of the atoms at various times
    for j in [10, 12, int(simulationLength/100)]:
        ax = plt.figure().add_subplot(projection='3d')
        ax.set_title('Position of atoms at time {}'.format(j * dt))
        for i in range(len(atoms)):
            name: str = 'Ar:{}'.format(i)
            ax.scatter(allAtoms[name]['position'][j][0], allAtoms[name]['position'][j][1], allAtoms[name]['position'][j][2])
        ax.set_xticks(np.arange(0, BOXSIZE, numberOfAtoms))
        ax.set_yticks(np.arange(0, BOXSIZE, numberOfAtoms))
        ax.set_zticks(np.arange(0, BOXSIZE, numberOfAtoms))
        ax.set_xlabel('x [Å]')
        ax.set_ylabel('y [Å]')
        ax.set_zlabel('z [Å]')
        plt.show()


    # Plot the position of 1 atom in the z-direction, located in the x-y plane.
    plt.title('Position of 1 atom in the z-direction')
    pos = np.array([allAtoms['Ar:{}'.format(0)]['position'][i][2] for i in range(len(allAtoms['Ar:0']['position']))])
    plt.plot(t, pos)
    plt.xlabel('Time [s]')
    plt.grid()
    plt.ylabel('Position [Å]')
    #plt.show()
    plt.savefig('pos-z_direction.png')
    plt.clf()

    # Plot the energy of the system
    plt.title('Energies of the system')
    plt.plot(t, np.array(Potential[1:]), label='$\mathcal{V}(t)$')
    K = 0
    for i in range(len(atoms)):
        vel = np.array([allAtoms['Ar:{}'.format(i)]['velocity'][j] for j in range(len(allAtoms['Ar:0']['position']))])
        K += 0.5 * atoms[i].mass * vel**2
    
    K = K[:,0] + K[:, 1] + K[:, 2]
    plt.plot(t, K , label='$\mathcal{K}(t)$')
    
    plt.plot(t, K + np.array(Potential[1:]), label='$\mathcal{H}(t)$')
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [eV]')
    plt.grid()
    #plt.show()
    plt.savefig('energy_gas.png')
    plt.clf()

    # Plot the temperature of the system
    plt.title('Temperature of the system')
    plt.plot(t, Temperature)
    plt.grid()
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.savefig('temperature.png')


if __name__ == '__main__':
    """
    state = State.bouncing
    runSimulation(state)
    """
    state = State.lattice
    runSimulation(state, 5)