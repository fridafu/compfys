from pylab import *

n = [10, 100, 1000]


for nn in n:

	f = zeros(nn)

	myfile = open('n%s.txt'%str(nn), 'r')
	i = 0
	for line in myfile:
		f[i] = float(line)
		i += 1


	x = linspace(0,1,nn)

	plot(x, f)
	hold('on')

xlabel('$x_i$')
ylabel('$v(x_i)$')
legend(['10', '100', '1000'])
show()
