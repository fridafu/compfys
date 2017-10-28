from pylab import *

myangle = open('relfile2.txt')
mytime = open('reltime2.txt')
myangle1 = open('relfile1.txt')
mytime1 = open('reltime1.txt')

t = []
r = []
t1 = []
r1 = []

for angle, time in zip(myangle, mytime):
	r.append(float(angle))
	t.append(float(time))

for angle, time in zip(myangle1, mytime1):
	r1.append(float(angle))
	t1.append(float(time))



myangle.close()
mytime.close()
myangle1.close()
mytime1.close()

t = array(t)
r = array(r)
t1 = array(t1)
r1 = array(r1)

times = linspace(0,100,10000)
m, c = polyfit(t,r,1)

times1 = times

m1, c1 = polyfit(t1, r1, 1)



angles = c + m*times

angles1 = c1 + m1*times

print(angles[-1])



print(angles1[-1])

plot(t1, r1)
plot(times, angles1, '--')
plot(t, r)
plot(times, angles, '--')
legend(['Perihelion angle with relaivistic correction', 'Linear fit of perihelion angle with relativistic correction', 'Perihelion angle with classical interaction', 'Linear fit of perihelion angle with classical interaction'])
xlabel('Time [years]')
ylabel('Angle [arc seconds]')
show()


