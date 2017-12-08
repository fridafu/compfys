from numpy import *
from pylab import *
# path = '/local/narga/os187/sprocess/'

files = ['test.dat']
# histogram for money distribution
dtype1 = [('counts', 'f8'), ('m','f8')]
# read in data
hist = loadtxt(files[0], dtype=dtype1)


y = hist['counts']
ddx = hist['m'][2] - hist['m'][1]


area = trapz(y, dx=ddx)

figure()
binno = size(hist['m'])


print(trapz(y, dx=ddx))

plot(hist['m'], hist['counts'])
#histogram(hist['m'], bins=binno)#hist['counts'])

show()
