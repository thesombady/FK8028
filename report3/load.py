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
        plt.savefig('energy{}.png'.format('bouncing'))
    else:
        plt.savefig('energy{}.png'.format('lattice'))
    #plt.show()
    plt.clf()


def plot_temperature(velocities):
    plt.title('Temperature of the system')
    dt = 1e-15
    ts = (np.array([i for i in range(len(velocities))]) * dt)[3000:]
    kinetic = np.array([sum(velocities[i]) for i in range(len(velocities))])
    print('Average kinetic energy: {}'.format(np.mean(kinetic)))
    temperature = kinetic * 2 / (3 * len(velocities[0]) * BOLTZMANN)
    temperature = temperature[3000:]
    average_temperature = np.mean(temperature)
    print('Average temperature: {}'.format(average_temperature))
    plt.plot(ts, temperature, label='$T$')
    plt.plot(ts, [average_temperature for _ in range(len(ts))], label='$\langle T \\rangle$', linestyle='--')
    plt.title('Temperature of the system')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    #plt.show()
    plt.savefig('temperature.png')
    plt.clf()

def plot_radial(radial):
    print(radial.shape)
    size = 18.09 # Å 
    dt = 1e-15
    rs = np.linspace(0, size, len(radial[0]))
    plt.title('Radial distribution function')
    for t in [2000, 5000, 10000, 15000]:
        plt.plot(rs, radial[t], label=f't={t * dt:.1e}')
    plt.legend()
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.grid()
    plt.savefig('radial.png')
    plt.clf()
    # Computing the average of the radial distribution function, after equilibration

    radial = radial[3000:]
    average_radial = np.mean(radial, axis=0)

    plt.title('Average radial distribution function')
    plt.plot(rs, average_radial)
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.grid()
    plt.savefig('average_radial.png')
    plt.clf()

    plt.title('Radial distribution function before equilibration')
    for t in [0, 200, 500, 1500]:
        plt.plot(rs, radial[t], label=f't={t * dt:.1e}')
    plt.legend()
    plt.ylabel('g(r)')
    plt.xlabel('r (Å)')
    plt.grid()
    plt.savefig('radial_before.png')
    plt.clf()

if __name__ == '__main__':
    #pos = load('lattice_pos.txt')
    vel = load_array('lattice_vel.txt')
    pot = load_array('lattice_pot.txt')
    potential = np.array([sum(pot[i]) for i in range(len(pot))])
    #plot_energy(pot, vel, 'Energy of the argon liquid')
    #plot_temperature(vel)
    #plot_energy(pot, vel)
    plot_energy(potential, vel, 'Energy of the argon liquid')
    radial = load_array('lattice_hist.txt')
    plot_radial(radial)
    """
        Compute the specific heat capacity, C_v, of the system.
        We skip the first 3000 steps.
    """
    vel = load_array('lattice_vel.txt')
    pot = load_array('lattice_pot.txt')
    potential = np.array([sum(pot[i]) for i in range(3000, len(pot))]) 
    plot_temperature(vel)
    kinetic = np.array([sum(vel[i]) for i in range(3000, len(vel))]) / 125 * 1.602 * 10 ** (-19)
    
    variance = np.var(kinetic)

    temperature = kinetic * 2 / (3 * BOLTZMANN * 1.602 * 10 ** (-19))

    average_temperature = np.mean(temperature)

    c_v = (2/ (3 * BOLTZMANN * 1.602 * 10 ** (-19)) - 4 * 125 * variance / (9 * (BOLTZMANN * 1.602 * 10 ** (-19)) ** 3 * average_temperature ** 2)) ** (-1)

    print('Specific heat capacity {}'.format(c_v * 6.022 * 10**(23)))

    print('Average potential energy: {}'.format(np.mean(potential[3000:])))


        

