from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as pl

x = []
y = []
z = []
with open('coordinates.txt') as file:
    next(file)
    for line in file:
        value = line.strip().split(' ')
        value = [float(i) for i in value]
        x.append(value[0])
        y.append(value[1])
        z.append(value[2])


fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('Distance [m]')
ax.set_ylabel('Distance [m]')
ax.set_zlabel('Distance [m]')

pl.show()
pl.clf()
