import os
import numpy as np
import matplotlib.pyplot as plt

def load_pos():
    path = 'bounce_pos.txt'
    with open(path, 'r') as f:
        content = f.read()

    pos_info = []
    content = content.split('\n')
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

def load_vel():
    path = 'bounce_vel.txt'
    with open(path, 'r') as f:
        content = f.read()

    pos_info = []
    content = content.split('\n')
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



def plot(pos):
    shape = pos.shape
    i = [i * 1e-15 for i in range(shape[0])]
    pos1 = [pos[i, 0, 0] for i in range(shape[0])]
    pos2 = [pos[i, 1, 0] for i in range(shape[0])]
    plt.plot(i, pos1, label='x1')
    plt.plot(i, pos2, label='x2')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    pos = load_pos()
    vel = load_vel()
    plot(vel)