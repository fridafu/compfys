from pylab import *

files = [open("L40.txt", "r"), open("L60.txt", "r"), open("L80.txt", "r"), open("L100.txt", "r")]
Ls = [40*40, 60*60, 80*80, 100*100]

E_array = []
E2_array = []
M_array = []
M2_array = []
M_abs = []

T = [];

for index, file in enumerate(files):
	lines=file.readlines()
	T = [];
	E = [];
	E2 = []; 
	M = [];
	M2 = [];
	Mabs = []
	for x in lines:
		T.append(float(x.split(' ')[0]))
		E.append(float(x.split(' ')[3]))
		E2.append(float(x.split(' ')[4]))
		M.append(float(x.split(' ')[5]))
		M2.append(float(x.split(' ')[6]))
		Mabs.append(float(x.split(' ')[7]))
	E = array(E)/Ls[index]
	E2 = array(E2)/Ls[index]
	M = array(M)/Ls[index]
	M2 = array(M2)/Ls[index]
	Mabs = array(Mabs)/Ls[index]

	E_array.append(E)
	E2_array.append(E2)
	M_array.append(M)
	M2_array.append(M2)
	M_abs.append(Mabs)

figure(1)
xlabel('Temperature T [kK/J]')
ylabel('<E> per spin [J]')
plot(T,E_array[0])
plot(T,E_array[1])
plot(T,E_array[2])
plot(T,E_array[3])
legend(['L=40','L=60', 'L=80','L=100'])
show()

figure(2)
xlabel('Magnetization M [unitless]')
ylabel('<E> per spin [J]')
plot(T,M_abs[0])
plot(T,M_abs[1])
plot(T,M_abs[2])
plot(T,M_abs[3])
legend(['L=40','L=60', 'L=80','L=100'])
show()

