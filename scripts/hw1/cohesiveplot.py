from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl

n = []
e = []
ep = []
with open('cohesive.txt') as file:
    next(file)
    for line in file:
        value = line.strip().split(' ')
        n.append(int(value[0]))
        e.append(float(value[1]))
        ep.append(float(value[2]))

a = 5.256*10**-10
x = [a*i**3 for i in n]

pl.plot(x, e, '.')
pl.plot(x, ep, '.')
pl.xlabel('Volume [m^3]')
pl.ylabel('Eth [eV]')
pl.grid(b=True, which='both')
pl.tight_layout()
pl.legend(['Non-periodic', 'Periodic'])
pl.show()
pl.clf()

print(x[-1])
print(e[-1])
print(ep[-1])
