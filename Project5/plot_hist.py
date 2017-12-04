from numpy import *
from pylab import *
# path = '/local/narga/os187/sprocess/'

files = ['hist.dat']
# histogram for money distribution
dtype1 = [('counts', 'f8'), ('m','f8')]
# read in data
hist = loadtxt(files[0], dtype=dtype1)


figure()
binno = size(hist['m'])
print binno
plot(hist['m'], hist['counts'])
#histogram(hist['m'], bins=binno)#hist['counts'])

show()
