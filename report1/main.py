import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
import matplotlib.pyplot as plt

@dataclass
class Atom: # Argon atom
    position: np.ndarray
    velocity: np.ndarray
    radius: float = 1.0 # Ångström
    mass: float = 39.948 * 931.394* 10**9 / (299_792_458)**2 # GeV/

def Normalize(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)

def LennardJones(pos1: np.ndarray, pos2: np.ndarray) -> (float, np.ndarray):
    """
       Computes the Lennards-Jones potential and force between two atoms. 
    """
    sigma: float = 3.40 
    epsilon: float = 0.0104
    r: np.ndarray = pos1 - pos2 # Relative position

    delta = (sigma / np.linalg.norm(r)) ** 6
    
    potential: float = 4 * epsilon * (delta ** 2 - delta)

    force: np.ndarray = 4 * epsilon * (12 * delta ** 2 - 6 * delta) * (r) / np.linalg.norm(r)**2

    return potential, force
    
def plotForcePotential(V: List[float], F_x: List[float], range_: np.ndarray) -> None:
    """
        Plots the force and potential given a potential and force in the X-direction.
        Also plots the force as a differece quotient.
    """
    plt.title('Interatomic force in the x direction')
    plt.plot(range_, F_x, label='${F_x}$')
    plt.xlabel('Distance [$\AA$]')
    plt.ylabel('Force [eV/$\AA$]')
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
    
    

def Initalize() -> List[Atom]:
    """
        Returns a list of two atoms seperated by 4 Ångström.
    """
    atoms: List[Atom] = [
        Atom(
            np.array([0,0,0]),
            np.array([0,0,0])
        ),
        Atom(
            np.array([4,0,0]),
            np.array([4,0,0])
        )
    ]

    return atoms



def VerletMethodPosition(atom: Atom, dt: float, force: np.ndarray):
    """
        Computes the new position of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit GeV/c^2
    """
    pos: np.ndarray = atom.position + atom.velocity * dt + force / (2 * atom.mass) * dt**2

    atom.position = pos


def VerletMethodVelocity(atom: Atom, dt: float, forceNew: np.ndarray, forceOld: np.ndarray):
    """
        Computes the new velocity of an atom using the Verlet method.
        Force is in unit eV/Å, dt is in unit s and mass is in unit GeV/c^2
    """
    vel: np.ndarray = atom.velocity + (forceNew/atom.mass  + forceOld/atom.mass) / (2) * dt

    atom.velocity = vel


def KineticEnergy(atom: Atom, velocity: np.ndarray) -> float:
    """
        Computes the kinetic energy of an atom.
    """
    return 0.5 * atom.mass * np.linalg.norm(velocity)**2


def computeHamiltonian(atoms: List[Atom], velocity: np.ndarray, potential: np.ndarray) -> float:
    """
        Computes the Hamiltonian of the system.
    """
    kinetic_part: float = 0
    for i in range(len(atoms)):
        kinetic_part += atoms[i].mass * np.linalg.norm(velocity[i])**2
    
    return kinetic_part + potential


def plotSimulation(t: np.ndarray, atom1_pos_x: List[float], atom2_pos_x: List[float],
                atom1_vel_x: List[float], atom2_vel_x: List[float], atom1_Kinetic: List[float],
                atom2_Kinetic: List[float], Potential: List[float], Hamiltonian: List[float]) -> None:

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
    plt.plot(t, Potential, label='V(t)')

    plt.plot(t, np.array(atom2_Kinetic) + np.array(atom1_Kinetic), label='$K_{tot}(t)$')

    plt.plot(t, Hamiltonian, label='$H(t)$')

    plt.grid()
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [eV]')
    plt.legend()
    plt.savefig('Different energies in the system.png')
    # plt.show()


def runSimulation(plotSim: bool, plotHamiltonian: bool, dt = 0.00001):
    """
        Run the simulation 
    """
    atoms = Initalize()

    T: float = 1
    t: float = 0

    """
    atoms_pos_x = {
        atoms[i]: atoms[i].position[0] for i in range(len(atoms))
    }
    """
    if plotHamiltonian:
        dt = [0.1, 0.05, 0.01, 0.005, 0.001]
        for i in range(len(dt)):
            Hamiltonian = runSimulation(False, False, dt[i])
            plt.plot(np.linspace(0, T, len(Hamiltonian)), Hamiltonian, label=f'dt = {dt[i]}')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [eV]')
        plt.title('Hamiltonian with different time-steps')
        plt.grid()

        plt.savefig('Hamiltonian with different time-steps.png')

        return
        
    
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
        Iterates over the atoms and computes the force and potential between them.
        It's only two atoms in the system, however it's written in a way such that
        it can compute the position of more items, it's just the plotting that varies.
    """
    while t < T:
        """
        for i in range(len(atoms)):
            
            potential: float = 0
            oldForce: np.ndarray = np.array([0,0,0])

            for j in range(len(atoms)):
                if i == j:
                    continue

                potential_, force_ = LennardJones(atoms[i].position, atoms[j].position)
                potential += potential_
                oldForce = oldForce + force_

            VerletMethodPosition(atoms[i], dt, oldForce)
            
        for i in range(len(atoms)):

            newForce: np.ndarray = np.array([0,0,0])

            for j in range(len(atoms)):
                if i == j:
                        continue

                _, force_ = LennardJones(atoms[i].position, atoms[j].position)
                newForce = newForce + force_

            VerletMethodVelocity(atoms[i], dt, oldForce, newForce)

            if i == 1:
                Potential.append(potential)
                atom1_pos_x.append(atoms[i].position[0])
                atom1_vel_x.append((atoms[i].velocity[0]))
                atom1_Kinetic.append(KineticEnergy(atoms[i], atoms[i].velocity[-1]))
            else:
                atom2_pos_x.append(atoms[i].position[0])
                atom2_vel_x.append((atoms[i].velocity[0]))
                atom2_Kinetic.append(KineticEnergy(atoms[i], atoms[i].velocity[-1]))
        
        Hamiltonian.append(computeHamiltonian(atoms, [atoms[i].velocity[-1], atoms[i].velocity[-1]], potential))
        print(f'Completed {t/T}%')
    
        # os.system('clear')
        """
        atom1 = atoms[0]
        atom2 = atoms[1]

        potential, force1 = LennardJones(atom1.position, atom2.position)
        force2 = -force1

        VerletMethodPosition(atom1, dt, force1)
        VerletMethodPosition(atom2, dt, force2)

        _, newForce = LennardJones(atom1.position, atom2.position)

        VerletMethodVelocity(atom1, dt, newForce, force1)
        VerletMethodVelocity(atom2, dt, newForce, force2)

        Potential.append(potential)

        atom1_pos_x.append(atom1.position[0])
        atom2_pos_x.append(atom2.position[0])
        atom1_vel_x.append((atom1.velocity[0]))
        atom2_vel_x.append((atom2.velocity[0]))

        atom1_Kinetic.append(KineticEnergy(atom1, atom1.velocity[-1]))
        atom2_Kinetic.append(KineticEnergy(atom2, atom2.velocity[-1]))

        Hamiltonian.append(computeHamiltonian(atoms, [atom1.velocity[-1], atom2.velocity[-1]], potential))


            
        t += dt

    t = np.linspace(0, T, len(atom1_pos_x)) # Time array

    if plotSim:
        plotSimulation(t, atom1_pos_x, atom2_pos_x, atom1_vel_x, atom2_vel_x, atom1_Kinetic, atom2_Kinetic, Potential, Hamiltonian)

    return Hamiltonian



if __name__ == '__main__':
    #computeForceAndPotentialForPlot()
    runSimulation(True, False)
