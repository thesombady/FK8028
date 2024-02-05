import numpy as np
from numpy import typing as npt
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt
from enum import Enum, auto

type Vector = npt.ArrayLike[npt.float64]

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

    spacing = BOXSIZE / (numberOfAtoms1D - 1) * 0.85

    print('Size: ' + str(BOXSIZE))

    counter = 0

    for i in range(numberOfAtoms1D):
        for j in range(numberOfAtoms1D):
            for k in range(numberOfAtoms1D):

                velocity = np.random.standard_normal(3) * 1e-10 # Å/s

                atoms.append(
                    Atom(
                        f'Ar:{counter}',
                        np.array([1 + (i) * spacing, 1 + (j) * spacing, 1 + (k) * spacing]),
                        velocity
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
    K = 0
    for v in velocities:
        K += 0.5 * np.linalg.norm(v)**2
    
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
        dt = 0.8e-15
    else:
        T = 1e-12
        dt = 1e-14
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

    Accelerations = []

    # Computes the initial acceleration for each atom to avoid recomputing it in the loop

    Potential: List[float] = []

    allAtoms = {atom.name: {'position': [], 'velocity': [], 'potential': []} for atom in atoms}

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
    
    while t < T:
        """
            Loop that iterates over the Velocity Verlet method
        """
        
        for i in range(len(atoms)):
            name: str = 'Ar:{}'.format(i)
            pos: Vector = VerletMethodPosition(atoms[i], dt, Accelerations[i])
            pos = validateBoundary(pos)
            atoms[i].position = pos 
            (allAtoms[name]['position']).append(atoms[i].position)
        
        newAccelerations = []

        computedValues = dict()

        for i in range(len(atoms)):
            name: str = 'Ar:{}'.format(i)
            
            potential: float = 0
            newForce: Vector = np.array([0.0, 0.0, 0.0])
            for j in range(len(atoms)):
                if i == j:
                    continue
                if computedValues.get((j, i)) is not None:
                    p, force_ = computedValues[(j, i)]
                    force_ *= -1
                else:
                    p, force_ = LennardJones(atoms[i].position, atoms[j].position)
                    computedValues[(i, j)] = (p, force_)
                potential += p
                newForce += force_
            
            Potential.append(potential)
            
            atoms[i].velocity = VerletMethodVelocity(atoms[i], dt, newForce/atoms[i].mass, Accelerations[i])

            (allAtoms[name]['velocity']).append(atoms[i].velocity)            

            newAccelerations.append(newForce / atoms[i].mass)

            (allAtoms[name]['potential']).append(potential)
        

        velo: List[Vector] = [allAtoms[name]['velocity'][-1] for name in allAtoms]

        K, Temp = kineticTemperature(velo)

        if t < 1e-13 and t % 50 == 0:

            for i in range(len(atoms)):
                velo[i] *= np.sqrt(94.4 / Temp)

            K, Temp = kineticTemperature(velo)

            for i in range(len(atoms)):
                atoms[i].velocity = velo[i]

        Accelerations = newAccelerations

        t += dt
        print(f'Completed {round(t/T, 4) * 100} %')

    t = np.linspace(0, T, len(allAtoms['Ar:0']['position']))

    if state == State.bouncing:
        """
            Plot the various plots for the Bouncing test.
        """
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
        plt.title('Energy of two particles')
        plt.plot(t, Kinetic, label="$\mathcal{K}(t)$")
        plt.plot(t, 1/2*Kinetic + potential, label="$\mathcal{H}(t)$")
        plt.plot(t, potential, label="$\mathcal{V}(t)$")
        plt.xlabel("Time [s]")
        plt.ylabel("Velocity [eV]")
        plt.grid()
        plt.legend()
        plt.savefig('bsE.png')
        #plt.show()
        plt.clf()
        return 
    
    simulationLength = len(allAtoms['Ar:0']['position'])
    print(T/dt)
    i = 0
    for j in [int(simulationLength * 1 / 8), int(simulationLength * 1 / 4), int(simulationLength * 1 / 2)]:

        ax = plt.figure().add_subplot(projection='3d')
        ax.set_title('Position of atoms at t = {} s'.format(j*dt))        
        ax.set_xticks(np.arange(0, BOXSIZE, numberOfAtoms))
        ax.set_yticks(np.arange(0, BOXSIZE, numberOfAtoms))
        ax.set_zticks(np.arange(0, BOXSIZE, numberOfAtoms))

        for i in range(len(atoms)):
            name = 'Ar:{}'.format(i)
            position = allAtoms[name]['position'][j]
            ax.scatter(position[0], position[1], position[2])
        
        ax.set_xlabel('x [Å]')
        ax.set_ylabel('y [Å]')
        ax.set_zlabel('z [Å]')
        #ax.set_clabel('z [Å]')

        #plt.savefig('lattice{}.png'.format(i))
        plt.show()
        i += 1
        

    plt.title('Position of 1 atom in the z-direction')
    pos = np.array([allAtoms['Ar:{}'.format(i)]['position'][3][2] for i in range(len(allAtoms['Ar:0']['position']))])
    plt.plot(t, pos)
    plt.xlabel('Time [s]')
    plt.ylabel('Position [Å]')
    plt.legend()
    #plt.show()
    plt.savefig('pos-z_direction.png')




if __name__ == '__main__':
    """
    state = State.bouncing
    runSimulation(state)
    """
    state = State.lattice
    #runSimulation(state, 5)


    potential, force = LennardJones(np.array([6.07,0.0,0.0]), np.array([0.07,0.0,0.0]))
    pos = VerletMethodPosition(Atom('Ar:0', np.array([0.07,0,0]), np.array([0,0,0])), 1e-15, force / ATOM_MASS)
    potential, newforce = LennardJones(pos, np.array([6.07,0.0,0.0]))
    vel = VerletMethodVelocity(Atom('Ar:0', np.array([0.07,0,0]), np.array([0,0,0])), 1e-15, newforce / ATOM_MASS, force / ATOM_MASS)
    #print(pos)
    #print(newforce)
    #print(vel)
    print(force)
    print(ATOM_MASS)