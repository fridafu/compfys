from numpy import *
from pylab import *
from math import gamma
# path = '/local/narga/os187/sprocess/'

files = ['Dalpha05_lambda0_N1000.dat', 'Dalpha1_lambda025_N1000.dat', 'Dalpha15_lambda025_N1000.dat', 'Dalpha2_lambda025_N1000.dat']
# histogram for money distribution
dtype1 = [('counts', 'f8'), ('m','f8')]
# read in data
X1 = loadtxt(files[0], dtype=dtype1)
X2 = loadtxt(files[1], dtype=dtype1)
X3 = loadtxt(files[2], dtype=dtype1)
X4 = loadtxt(files[3], dtype=dtype1)

step = 100
M1 = X1['m'][1::step]
C1 = array([sum(X1['counts'][i:i+step]) for i in range(0, len(X1['counts']), step)])

M2 = X2['m'][1::step]
C2 = array([sum(X2['counts'][i:i+step]) for i in range(0, len(X2['counts']), step)])
M3 = X3['m'][1::step]
C3 = array([sum(X3['counts'][i:i+step]) for i in range(0, len(X3['counts']), step)])
M4 = X4['m'][1::step]
C4 = array([sum(X4['counts'][i:i+step]) for i in range(0, len(X4['counts']), step)])



norm = sum(C1)
figure()

loglog(M1, C1/norm, 'r')
hold('on')

"""
loglog(M2,C2/norm, 'g')
hold('on')

loglog(M3,C3/norm, 'b')
hold('on')

loglog(M4,C4/norm, 'k')
hold('on')
"""




w = 10**12*M1**(-(8))




loglog(M1, w)




"""
legend([r'$\gamma$ = 0',r'$\gamma$ = 1',r'$\gamma$ = 2',r'$\gamma$ = 3',r'$\gamma$ = 4', r'Power Law'])
"""

legend([r'$\gamma$ = 1', r'Power Law'])

xlabel('Money [bitcoins]')
ylabel('Probability')
show()


#x = X1['m']/2.0
#omega = x**(-0.75)
#n = 1 + 3*0.25/(0.75)
#an = (n**n/(math.gamma(n)))*0.026
#Pnx = an*x**(n-1)*exp(-n*x)
#f = 0.5*exp(-0.5*x)
"""
start = 1
step = 1
norm = sum(X1['counts'])
figure()
loglog(X1['m'][start::step], X1['counts'][start::step]/norm, 'r')
hold('on')





loglog(X2['m'][start::step], X2['counts'][start::step]/norm, 'g')
hold('on')
loglog(X3['m'][start::step], X3['counts'][start::step]/norm, 'b')
hold('on')
loglog(X4['m'][start::step], X4['counts'][start::step]/norm, 'k')
hold('on')
loglog(X5['m'][start::step], X5['counts'][start::step]/norm, 'm')
hold('on')

w = X1['m']**(-(1 + 1))


loglog(X1['m'], w)
hold('on')
w = X1['m']**(-(1 + 2))
loglog(X1['m'][start::step], w[start::step])




#hold('on')
#plot(lam025['m'], Pnx, 'k')
xlabel('Money [bitcoins]')
ylabel('Probability')
legend([r'$\gamma$ = 0',r'$\gamma$ = 1',r'$\gamma$ = 2',r'$\gamma$ = 3',r'$\gamma$ = 4'])
#xlim([0,10])
#histogram(hist['m'], bins=binno)#hist['counts'])

show()
"""

