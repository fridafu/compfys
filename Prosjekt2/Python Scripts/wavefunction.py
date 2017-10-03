# 2e) plotting the wavefunction

import matplotlib.pyplot as plt
import numpy as np


filenames = ['wr001.txt','wr05.txt','wr1.txt','wr5.txt']
omega = [0.01, 0.5, 1, 5]
rho_max = [50, 8, 5, 2.5]


def read_file(filename, n=100):
    E =np.zeros((n,3))
    for i, line in enumerate(open(filename, 'r')): 
        line =line.rstrip() # remove whitespaces
        if line:            # include only non-empty lines.
            ref = line.split();
            E[:,i] = [float(x) for x in ref]
        E0 = E[:,0]; E1 = E[:,1];E2 = E[:,2]
    return E0, E1, E2

def plotting(E0, E1, E2, rhomax,wr,interact = 1, n=100):
    h = rhomax/float(n)
    rho = np.linspace(h,rhomax,n)
    #plot the normalized wave functions for different values of wr
    plt.plot(rho,(E0**2)/np.linalg.norm(E0**2),rho,(E1**2)/np.linalg.norm(E1**2), rho, (E2**2)/np.linalg.norm(E2**2))
    plt.xlabel(r'$\rho$',fontsize = 30), plt.ylabel(r'$|\psi(\rho)|^2$',fontsize = 30)
    plt.axis(fontsize=30)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    #change title if there is interaction
    if interact == 1:
        plt.title(r'$Interaction,$ $\omega_r = $'+str(wr)+ r'$,\,\rho_{max} = $'+str(rhomax),fontsize = 30)
    else:
         plt.title(r'$No\,interaction,$ ' r'$\rho_{max} = $'+str(rhomax),fontsize = 30)
    plt.legend([r'$E_0$',r'$E_1$',r'$E_2$'], fontsize = 30)
    plt.show()

for i in range(len(filenames)):
    E0, E1, E2 = read_file(filenames[i])
    plotting(E0, E1, E2,rho_max[i], omega[i])

def waves(E0C, E0noC, rhomax,n=100):# plot the different wavefunctions E0 w/ and wo/ the interaction
    h = rhomax/float(n)
    rho = np.linspace(h,101*h,n)
    plt.plot(rho,(E0C**2)/np.linalg.norm(E0C**2), rho,(E0noC**2)/np.linalg.norm(E0noC**2))
    plt.xlabel(r'$\rho$',fontsize = 30), plt.ylabel(r'$|\psi(\rho)|^2$',fontsize = 30)
    plt.axis(fontsize=30)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.title(r'With and without $\frac{1}{\rho}$ $\omega_r = 1$ $\rho_{max} = 5$',fontsize = 30)
    plt.legend([r'$E_0$ with interaction','$E_0$ without interaction'], fontsize = 30)
    plt.show()


E0C, E1x, E2x = read_file('wr1.txt')

E0noC, E1nx, E2nx = read_file('wr1nointeract.txt')
plotting(E0noC, E1nx, E2nx, 5, 1, interaction = 0)
waves(E0C, E0noC, rhomax=5)
