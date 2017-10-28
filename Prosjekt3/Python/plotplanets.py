from pylab import *

mypositionfile = open("solarsystem.txt", "r")


r = []
u = 0

for linep in mypositionfile:
	if u%10 == 0:

		numplanets = int(len(linep.split())/3)
		start = 0
		end = 3
		templist = []

		for i in range(numplanets):
			templist.append(list(map(float, linep.split()))[start:end])
			start += 3
			end += 3
	
		r.append(templist)
	u += 1
	
mypositionfile.close()


planets =['Earth', 'Mars', 'Saturn', 'Uranus', 'Jupiter', 'Venus', 'Mercury', 'Neptune', 'Pluto', 'Sun']




r = array(r) #(timestep, planet, coordinate)


for i in range(len(r[0,:,0])):
	plot(r[:,i,0], r[:,i,1])
	
lim = 70

ylim([-lim, lim])
xlim([-lim,lim])


legend(planets)
xlabel('x [AU]')
ylabel('y [AU]')

show()
