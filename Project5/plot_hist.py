from numpy import *
from pylab import *
# path = '/local/narga/os187/sprocess/'

files = ['test.dat']
# histogram for money distribution
dtype1 = [('counts', 'f8'), ('m','f8')]
# read in data
hist = loadtxt(files[0], dtype=dtype1)
print(hist['m'])
print(hist['counts'])


figure()
binno = size(hist['m'])
print(binno)
plot(hist['m'], hist['counts'])
#histogram(hist['m'], bins=binno)#hist['counts'])

show()
