import os
import numpy as np
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6 # eV/c^2

ATOM_MASS = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2


def load(path: str):
    with open(path, 'r') as f:
        content = f.read()

    pos_info = []
    content = content.split('\n')[:-1]
    for i in range(len(content)):
        element = content[i].split('\t')
        rows = []
        for j in range(len(element)):
            cellis = []
            cell = (element[j].replace('[', '').replace(']', '')).split(',')
            for i in range(3):
                cellis.append(float(cell[i]))
            rows.append(cellis)
        pos_info.append(np.array(rows))

    
    return np.array(pos_info)

def load(path: str):
    with open(path, 'r') as f:
        content = f.read()

    pos_info = []
    content = content.split('\n')[:-1]
    for i in range(len(content)):
        element = content[i].split('\t')
        rows = []
        for j in range(len(element)):
            cellis = []
            cell = (element[j].replace('[', '').replace(']', '')).split(',')
            for i in range(3):
                cellis.append(float(cell[i]))
            rows.append(cellis)
        pos_info.append(np.array(rows))

    
    return np.array(pos_info)

def load_single(path: str):
    with open(path, 'r') as f:
        content = f.read()

    pos_info = []
    content = content.split('\n')[:-1]
    potential = []
    for i in range(len(content)):
        element = content[i].split('\t')
        total_potential = 0
        for j in range(len(element)):
            total_potential += float(element[j])

        potential.append(total_potential)
    
    return np.array(potential)

def load_temp(path: str):
    with open(path, 'r') as f:
        content = f.read()

    content = content.split('\n')[:-1]

    return np.array([float(i) for i in content])

def plot(pos):
    shape = pos.shape
    i = [i * 1e-15 for i in range(shape[0])]
    pos1 = [pos[i, 0, 0] for i in range(shape[0])]
    pos2 = [pos[i, 1, 0] for i in range(shape[0])]
    plt.plot(i, pos1, label='x1')
    plt.plot(i, pos2, label='x2')
    plt.legend()
    plt.show()


def plot_energy(potential, velocity):
    plt.title('Energy of the system')
    ts = np.array([i for i in range(len(potential))]) * 1e-15
    plt.plot(ts, potential, label='$\mathcal{V}$')
    kinetic = []
    for i in range(len(velocity)):
        vel = 0.0
        for j in range(len(velocity[i])):
            vel += np.linalg.norm(velocity[i][j])**2
        kinetic.append(vel)
    kinetic = np.array(kinetic) * 0.5 * ATOM_MASS
    plt.plot(ts, kinetic, label='$\mathcal{K}$')
    plt.plot(ts, 2 * kinetic + potential, label='$\mathcal{H}$')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.show()


def plot_temperature(temperature):
    plt.title('Temperature of the system')
    ts = np.array([i for i in range(len(temperature))]) * 1e-15
    plt.plot(ts[2000:], temperature[2000:], label='Temperature')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.savefig('temperature.png')
    #plt.show()


if __name__ == '__main__':
    #pos = load('lattice_pos.txt')
    #vel = load('lattice_vel.txt')
    #pot = load_single('lattice_pot.txt')
    temp = load_temp('lattice_temp.txt')
    #plot_energy(pot, vel)
    plot_temperature(temp)