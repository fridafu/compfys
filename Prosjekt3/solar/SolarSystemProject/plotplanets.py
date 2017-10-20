from pylab import *

mypositionfile = open("solarsystem.txt", "r")
mytimefile = open("times.txt", "r")

r = []
t = []

for linep, linet in zip(mypositionfile, mytimefile):
	numplanets = int(len(linep.split())/3)
	start = 0
	end = 3
	templist = []
	for i in range(numplanets):
		templist.append(list(map(float, linep.split()))[start:end])
		start += 3
		end += 3
	
	r.append(templist)
	t.append(float(linet))



planets =['Earth', 'Mars', 'Saturn', 'Uranus', 'Jupiter', 'Venus', 'Mercury', 'Neptune', 'Pluto', 'Sun']
mypositionfile.close()
mytimefile.close()

r = array(r) #(timestep, planet, coordinate)
print(r.shape)
t = array(t)


for i in range(len(r[0,:,0])):
	plot(r[:,i,0], r[:,i,1])
	hold('on')
lim = 40

ylim([-lim, lim])
xlim([-lim,lim])


legend(planets)
xlabel('x [AU]')
ylabel('y [AU]')

show()
