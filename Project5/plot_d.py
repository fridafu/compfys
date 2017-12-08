from numpy import *
from pylab import *
from math import gamma
# path = '/local/narga/os187/sprocess/'

files = ['Dalpha1_lambda025_N1000.dat', 'Ealpha1_lambda025_gamma_1_N1000.dat','Ealpha1_lambda025_gamma_2_N1000.dat','Ealpha1_lambda025_gamma_3_N1000.dat','Ealpha1_lambda025_gamma_4_N1000.dat']
# histogram for money distribution
dtype1 = [('counts', 'f8'), ('m','f8')]
# read in data
X1 = loadtxt(files[0], dtype=dtype1)
X2 = loadtxt(files[1], dtype=dtype1)
X3 = loadtxt(files[2], dtype=dtype1)
X4 = loadtxt(files[3], dtype=dtype1)
X5 = loadtxt(files[4], dtype=dtype1)

#x = X1['m']/2.0
#omega = x**(-0.75)
#n = 1 + 3*0.25/(0.75)
#an = (n**n/(math.gamma(n)))*0.026
#Pnx = an*x**(n-1)*exp(-n*x)
#f = 0.5*exp(-0.5*x)

norm = sum(X1['counts'])
figure()
loglog(X1['m'], X1['counts']/norm, 'r')
hold('on')
loglog(X2['m'], X2['counts']/norm, 'g')
hold('on')
loglog(X3['m'], X3['counts']/norm, 'b')
hold('on')
loglog(X4['m'], X4['counts']/norm, 'k')
hold('on')
loglog(X5['m'], X5['counts']/norm, 'm')

#hold('on')
#plot(lam025['m'], Pnx, 'k')
xlabel('Money [bitcoins]')
ylabel('Probability')
legend([r'$\gamma$ = 0',r'$\gamma$ = 1',r'$\gamma$ = 2',r'$\gamma$ = 3',r'$\gamma$ = 4'])
#xlim([0,10])
#histogram(hist['m'], bins=binno)#hist['counts'])

show()
