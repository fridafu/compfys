from pylab import *

files = [open("L40.txt", "r"), open("L60.txt", "r"), open("L80.txt", "r"), open("L100.txt", "r")]
Ls = [40*40, 60*60, 80*80, 100*100]

k = 1.38064852e-23

E_array = []
E2_array = []
M_array = []
M2_array = []
M_abs = []
cv_array = []
sus_array = []
T = []

for i, file in enumerate(files):
    lines=file.readlines()
    T = []
    E = []
    E2 = []
    M = []
    M2 = []
    Mabs = []
    cv = []
    sus = []

    for x in lines:
        T.append(float(x.split(' ')[0]))
        E.append(float(x.split(' ')[3]))
        E2.append(float(x.split(' ')[4]))
        M.append(float(x.split(' ')[5]))
        M2.append(float(x.split(' ')[6]))
        Mabs.append(float(x.split(' ')[7]))
#print E2[-1] - (E[-1])**2)/(k*(T[-1])**2)
        cv.append((E2[-1] - (E[-1])**2)/((T[-1])**2))
        sus.append((M2[-1] - (Mabs[-1])**2)/(T[-1]))

    E = array(E)/Ls[i]
    E2 = array(E2)/Ls[i]
    M = array(M)/Ls[i]
    M2 = array(M2)/Ls[i]
    Mabs = array(Mabs)/Ls[i]
    cv = array(cv)/Ls[i]
    sus = array(sus)/Ls[i]
    
    E_array.append(E)
    E2_array.append(E2)
    M_array.append(M)
    M2_array.append(M2)
    M_abs.append(Mabs)
    cv_array.append(cv)
    sus_array.append(sus)

figure(1)
xlabel('Temperature T [kK/J]')
ylabel('<E> per spin [J]')
plot(T,E_array[0])
plot(T,E_array[1])
plot(T,E_array[2])
plot(T,E_array[3])
legend(['L=40','L=60', 'L=80','L=100'])

figure(2)
xlabel('Temperature T[kK/J]')
ylabel('Magnetization <|M|> per spin [unitless]')
plot(T,M_abs[0])
plot(T,M_abs[1])
plot(T,M_abs[2])
plot(T,M_abs[3])
legend(['L=40','L=60', 'L=80','L=100'])

figure(3)
xlabel('Temperature T [kK/J]')
ylabel('Heat capacity $<C_v>$ per spin')
plot(T,cv_array[0])
plot(T,cv_array[1])
plot(T,cv_array[2])
plot(T,cv_array[3])
legend(['L=40','L=60', 'L=80','L=100'])

figure(4)
xlabel('Temperature T [kK/J]')
ylabel('Magnetic Susceptibility $<\chi>$ per spin')
plot(T,sus_array[0])
plot(T,sus_array[1])
plot(T,sus_array[2])
plot(T,sus_array[3])
legend(['L=40','L=60', 'L=80','L=100'])
show()

A = matrix([[1, 1/100],[1, 1/80],[1,1/60],[1,1/40]])
T = matrix([[2.275],[2.280],[2.285],[2.290]])

x = A.getI() * T
print(x)

