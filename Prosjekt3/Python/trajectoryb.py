# plotting the trajectory for Sun and Earth
# forward Euler
import numpy as np
import matplotlib.pyplot as plt

# reading in the coordinates from file

def read_xyz(filename,n):
    xy = np.zeros((n,2))
    for i, line in enumerate(open(filename, 'r')): 
        line =line.rstrip() # remove whitespaces
        if line:            # include only non-empty lines.
            ref = line.split()
            xy[i,:] = [float(j) for j in ref]
    return xy

def plot_trajectory(x, y):
    plt.hold('on')
    plt.scatter(0,0,c='y')
    plt.plot(x,y)
    plt.axis('equal')
    

xy = read_xyz('fetest1.txt', n=300000)
plot_trajectory(xy[:,0], xy[:,1])
xy_vv = read_xyz('vvtest1.txt', n=300000)
plot_trajectory(xy_vv[:,0], xy_vv[:,1])
#plt.legend['fe', 'vv']
plt.show()
