import numpy as np
import matplotlib.pyplot as plt

import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6 # eV/c^2

ATOM_MASS = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2

BOLTZMANN = 8.617333262145 * 1e-5 # eV/K
TEMPERATURE = 94.4 # K
BETA = 1 / (BOLTZMANN * TEMPERATURE) # 1 / eV


def load(path: str) -> np.ndarray:
    """
        Load the data from a file and return it as a numpy array
        : path { str } - The path to the file
        ; return { np.array } - The data from the file as a vector
    """
    with open(path, 'r') as file:
        content = file.read()
    content = content.split('\n')[:-1] 

    elements = []

    for i in range(len(content)):
        elements.append(float(content[i]))

    return np.array(elements) 


def load_array(path: str) -> np.ndarray:
    """
        Load the data from a file and return it as a numpy array
        : path { str } - The path to the file
        ; return { np.array } - The data from the file as a vector
    """
    with open(path, 'r') as file:
        content = file.read()
    content = content.split('\n')[:-1]
    
    out = np.zeros((len(content), 100))

    for i in range(len(content)):
        row = content[i].split('\t')
        out_row = np.zeros((100))
        
        for j in range(len(row)):
            out_row[j] = float(row[j])
        out[i] = out_row
    
    return out
    

def plot_radial(radial: np.ndarray):
    """
        Computes the radial function and plots it.
        : radial { np.ndarray } - The radial distribution function
    """
    
    rs = np.linspace(0, 18.04, 100) # Bin size

    average_radial = np.mean(radial, axis=0)

    plt.title('Average radial distribution function')
    #plt.plot(rs, radial[0], label = '$g_0(r)$')
    plt.plot(rs, average_radial, label = '$ \langle g(r)\\rangle = \\frac{1}{N}\\sum_{i = 1}^Ng_i(r)$')
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.grid()
    plt.legend()
    plt.savefig('average_radial.png')
    plt.clf()

def return_block(data: np.array, k: int):
    """
        Returns a block of data from the original data array.
        : data { np. array } - The data of which we want to compute the block average
        : k { int } - The size of the block, i.e. block-size
    """

    nb = int(np.floor(len(data) / k))

    average = np.mean(data)

    blocks = []
    
    for i in range(0,len(data), k):
        value = 1/k * sum(data[i:i+k])
        blocks.append(value)
    
    variance_block = 0.0

    block_average = np.mean(blocks)

    for i in range(nb):
        variance_block += 1/nb * (blocks[i] - average) ** 2
    
    std_err = np.sqrt(variance_block / (nb - 1))

    err = std_err / np.sqrt(2 * (nb - 1))
    return std_err, err, block_average

def make_figure(data: np.array, title: str, save_str: str):
    """
        Make a figure of the standard error and the mean value of the data.
        : data { np.array } - The data of which we want to compute the block average
        : title { str } - The title of the plot
        : save_str { str } - The name of the file to save the plot
    """
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    end = 16
    l1 = ax[1].plot(np.arange(1, end), np.ones(end - 1) * np.mean(data), color='blue', label='$\\bar{A}$')
    for i in range(1, end):
        k = 2 ** i 
        std_err, err, block_average = return_block(data, k)
        ax[0].errorbar(k, std_err, yerr = err, color='black', markersize = 2, fmt='o')
        
        l2 = ax[1].errorbar(k, block_average, yerr = std_err, color = 'black', markersize = 3, fmt='o', label='Block average' if i == 1 else '')

        if i == 15:
            print()
            print('sigma(V_k) = {}'.format(std_err))

            #print('Standard error {}\nand relative error {} of {}\n Mean {}'.format(std_err, err, title, block_average))
    
    ax[1].legend()

    ax[0].set_xlabel('Block length ')
    ax[0].set_ylabel('Standard deviation $\\sigma(A_k)$') 
    ax[0].set_title('(a) Standard deviation of {}'.format(title))
    ax[1].set_title('(b) Mean value of of {}'.format(title))
    ax[1].set_xlabel('Block length')
    ax[1].set_ylabel('Mean value $\\bar{A}_k$') 
 
    plt.savefig('{}.png'.format(save_str))
    #plt.show()
    plt.clf()



def plot_potential(potential: np.ndarray):
    """
        Plot the potential energy of the system
        : potential { np.ndarray } - The potential energy of the system
    """

    nu_var = np.var(potential[100_000:] / 125)


    mean_potential = np.mean(potential)

    c_v = 3 * BOLTZMANN / 2 + 125 * nu_var / ( BOLTZMANN * TEMPERATURE ** 2)
    print('Specific heat-capacity: {} J / (Mol K)'.format(c_v * 6.022 * 1e23 * 1.602 * 1e-19)) # J / (Mol K)


    print('Mean potentail: {} eV'.format(mean_potential))

    var_potential = np.var(potential)

    print('Error potentail: {} '.format(var_potential))

    print('Standard divation {} eV'.format(np.sqrt(var_potential)))

    iterations = np.arange(0, len(potential))

    plt.plot(iterations, potential, label = '$\mathcal{V}$')
    plt.plot(iterations, mean_potential * np.ones(len(potential)), label = '$\\bar{\mathcal{V}}$')
    plt.vlines(100_000, -15.5, -12.5, color = 'black', label = 'Equilibration')
    plt.title('Potential energy of the system')
    plt.xlabel('Number of iterations')
    plt.ylabel('$\mathcal{V}$ eV')
    plt.legend()
    plt.savefig('potential_energy.png')
    plt.clf()
    


if __name__ == '__main__':
    g_r = load_array('lattice_hist.txt')
    potential = load('lattice_pot.txt')
    plot_radial(g_r)
    make_figure(potential, '$A = \mathcal{V}$', 'var_p')
    plot_potential(potential)