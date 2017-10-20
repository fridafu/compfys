from pylab import *


myfilepos = open("solarsystem.txt", 'r')
myfilet = open("times.txt")


r = []
t = []

for linep, linet in zip(myfilepos, myfilet):
	r.append(list(map(float, linep.split())))
	t.append(float(linet))

r = array(r)
t = array(t)

myfilepos.close()
myfilet.close()

plot(r[:,0], r[:,1])
show()



	
	