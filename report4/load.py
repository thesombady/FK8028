import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List

SPEED_OF_LIGHT = 299_792_458 # m/s
GEV_CONVERSION = 931.394 * 10**6 # eV/c^2

ATOM_MASS = 39.948 * GEV_CONVERSION / (SPEED_OF_LIGHT * 1e10)**2 # eVs^2/Å^2

BOLTZMANN = 8.617333262145 * 10**-5 # eV/K



def load(path: str) -> np.ndarray:
    """
        Load the data from a file and return it as a numpy array
        : pat { str } - The path to the file
        ; return { np.array } - The data from the file as a vector
    """
    with open(path, 'r') as file:
        content = file.read()
    content = content.split('\n')[:-1] 

    elements = []

    for i in range(len(content)):
        elements.append(float(content[i]))
    return np.array(elements) 


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
    l1 = ax[1].plot(np.arange(1, 14), np.ones(13) * np.mean(data), color='blue', label='$\\bar{A}$')
    for i in range(1, 13):
        k = 2 ** i 
        std_err, err, block_average = return_block(data, k)
        ax[0].errorbar(i, std_err, yerr = err, color='black', markersize = 2, fmt='o')
        
        l2 = ax[1].errorbar(i, block_average, yerr = std_err, color = 'black', markersize = 3, fmt='o', label='Block average' if i == 1 else '')

        if i == 6:
            print('Standard error {} and relative error {} of {}, Mean {}'.format(std_err, err, title, block_average))
    
    ax[1].legend()

    ax[0].set_xlabel('Block size $(2^i)$')
    ax[0].set_ylabel('Standard deviation $\\sigma(A_k)$') 
    ax[0].set_title('(a) Standard deviation of {}'.format(title))
    ax[1].set_title('(b) Mean value of of {}'.format(title))
    ax[1].set_xlabel('Block size $(2^i)$')
    ax[1].set_ylabel('Mean value $\\bar{A}_k$') 
 
    plt.savefig('{}.png'.format(save_str))
    #plt.show()
    plt.clf()



if __name__ == '__main__':
    """
        N : Number of atoms, which we know
        Compute the block average for all the observables in the system
        We remove the first 3000 time-steps to remove the velocity rescaling and allow for the system to equilibrate
    """
    N = 125
    kinetic = load('lattice_kin.txt')
    potential = load('lattice_pot.txt')
    temperature = kinetic[3000:] * 2 / (3 * N * BOLTZMANN)
    make_figure(kinetic[3000:], '$A = \mathcal{K}$', 'var_k')
    make_figure(potential[3000:], '$A = \mathcal{V}$', 'var_p',)
    make_figure(temperature, '$A = T$', 'var_t')
    print('Potential average: {}'.format(np.mean(potential[3000:])))
    print('Kientic average: {}'.format(np.mean(kinetic[3000:])))
    print('Potential average: {}'.format(np.mean(temperature)))

