
import matplotlib
matplotlib.use('AGG')
from numpy import *
from matplotlib.pyplot import *
filenames = ['vescmin500hd.txt'] #,'vescmin500.txt','vescmin500_001.txt']#,'vescmin500_0001.txt']#,'vescplus50.txt','vescmin50.txt']# ['vesc100.txt','vescplus100.txt', 'vescmin100.txt']#['vesc.txt', 'vescplus.txt','vescmin.txt']#, 'beta25.txt','beta29.txt','beta299.txt','beta3.txt']#['vescplus1000.txt', 'vescminus1000.txt', 'vesc1000.txt']

def read_file(filename):
    
    mypositionfile = open(filename, 'r')
    mytimefile = open('times500hd.txt', 'r')
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

    mypositionfile.close()
    mytimefile.close()
    r = array(r) #(timestep, planet, coordinate)
    
    t = array(t)
    x = r[:,0,0]; y = r[:,0,1]
    
    return x, y

fig=figure() 
    
xlabel(r'x [AU]')
ylabel(r'y [AU]')

for i in range(len(filenames)):
    x, y = read_file(filenames[i])
    plot(x,y)
  
scatter(0,0,c='y')
axis('equal')
xlim((-45,1.5))
ylim((-10,20))
legend([r'$v_y = 2\sqrt{2}\pi$',r'$v_y = 2\sqrt{2}\pi-0.1$',r'$v_y = 2\sqrt{2}\pi-0.01$', r'Sun'])#,r'$\beta = 2.9$',r'$\beta =2.999$',r'$\beta = 3$',r'Sun'],loc='upper left')
savefig('escape500hdyears')
