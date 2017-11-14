from pylab import *
filen = 'L40.txt'
T = loadtxt(filen, usecols = 0)
L = loadtxt(filen, usecols = 1)[0]

E = loadtxt(filen,usecols = 3)

E2 = loadtxt(filen,usecols = 4)

M = loadtxt(filen,usecols = 5)

M2 = loadtxt(filen, usecols = 6)

absM = loadtxt(filen, usecols = 7)

figure(1)
subplot(211)

plot(T, E2 - E*E)

subplot(212)
plot(T, M2 - absM*absM)

show() 
