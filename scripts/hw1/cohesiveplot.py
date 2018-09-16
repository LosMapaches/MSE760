from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl

n = []
e = []
with open('cohesive.txt') as file:
    for line in file:
        value = line.strip().split(' ')
        n.append(int(value[0]))
        e.append(float(value[1]))

a = 5.256*10**-10
x = [a*i**3 for i in n]

pl.plot(x, e, '.')
pl.xlabel('Volume [m^3]')
pl.ylabel('Eth [eV]')
pl.grid(b=True, which='both')
pl.tight_layout()
pl.show()

print(x[-1])
print(e[-1])
