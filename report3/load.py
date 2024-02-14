import os
import numpy as np
import matplotlib.pyplot as plt

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6 # eV/c^2

ATOM_MASS = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2

BOLTZMANN = 8.617333262145 * 10**-5 # eV/K


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

def load_array(path: str):
    with open(path, 'r') as f:
        content = f.read()

    content = content.split('\n')[:-1]
    rows = []
    for i in range(len(content)):
        element = content[i].split('\t')
        row = []
        for j in range(len(element)):
            row.append(float(element[j]))
        rows.append(np.array(row))
    
    return np.array(rows)


def plot(pos):
    shape = pos.shape
    i = [i * 1e-15 for i in range(shape[0])]
    pos1 = [pos[i, 0, 0] for i in range(shape[0])]
    pos2 = [pos[i, 1, 0] for i in range(shape[0])]
    plt.plot(i, pos1, label='x1')
    plt.plot(i, pos2, label='x2')
    plt.legend()
    plt.show()


def plot_energy(potential, velocity, title):
    dt = 1e-15
    plt.title(title)
    ts = np.array([i for i in range(len(potential))]) * dt
    kinetic = np.array([sum(velocity[i]) for i in range(len(velocity))])
    plt.plot(ts, potential, label='$\mathcal{V}$')
    plt.plot(ts, kinetic, label='$\mathcal{K}$')
    plt.plot(ts, 2 * kinetic + potential, label='$\mathcal{H}$')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    if 'bouncing' in title:
        #plt.savefig('energy{}.png'.format{bouncing})
        ...
    else:
        #plt.savefig('energy{}.png'.format{lattice})
        ...
    plt.show()


def plot_temperature(velocities):
    plt.title('Temperature of the system')
    dt = 1e-15
    ts = np.array([i for i in range(len(velocities))]) * dt
    kinetic = np.array([sum(velocities[i]) for i in range(len(velocities))])
    temperature = kinetic * 2 / (3 * len(velocities[0]) * BOLTZMANN)
    plt.plot(ts[2000:], temperature[2000:], label='$T$')
    plt.title('Temperature of the system')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    #plt.show()
    plt.savefig('temperature.png')

def plot_radial(radial):
    print(radial.shape)
    size = 18.09 # Å 
    dt = 1e-15
    rs = np.linspace(0, size, len(radial[0]))
    plt.title('Radial distribution function')
    for t in [0, 2000, 5000, 10000, 15000]:
        plt.plot(rs, radial[t], label=f't={t * dt:.1e}')
    plt.legend()
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.grid()
    plt.savefig('radial.png')

if __name__ == '__main__':
    #pos = load('lattice_pos.txt')
    vel = load_array('bouncing_vel.txt')
    pot = load_array('bouncing_pot.txt')
    pot = np.array([sum(pot[i]) for i in range(len(pot))])
    plot_energy(pot, vel, 'Energy of the bouncing test')
    #plot_temperature(vel)
    #temp = load_temp('lattice_temp.txt')
    #plot_energy(pot, vel)
    #plot_temperature(temp)
    #radial = load_array('lattice_hist.txt')
    #plot_radial(radial)