# 2e) plotting the wavefunction

import matplotlib.pyplot as plt
import numpy as np

n = 100
#filename = raw_input('enter filename ')
E =np.zeros((n,3))

for i, line in enumerate(open('wr001.txt', 'r')): 
    line =line.rstrip() # remove whitespaces
    if line:            # include only non-empty lines.
        ref = line.split();
        E[:,i] = [float(x) for x in ref]
        
E0 = E[:,0]; E1 = E[:,1];E2 = E[:,2]
rhomax = float(raw_input('rho_max = '))
h = rhomax/n
rho = np.linspace(h,101*h,n)

#titles = raw_input('plot title')
plt.plot(rho,E0**2,rho,E1**2, rho, E2**2)
plt.xlabel(r'$\rho$',fontsize = 30), plt.ylabel(r'$|\psi(\rho)|^2$',fontsize = 30)
plt.axis(fontsize=30)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.legend([r'$E_0$',r'$E_1$',r'$E_2$'], fontsize = 30)
plt.title(r'Interaction, $\omega_r = 0.01,\,\rho_{max} = 50$',fontsize = 35)
plt.show()


