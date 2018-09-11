from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

x = []
y = []
z = []
with open('coordinates.txt') as file:
    for line in file:
        value = line.strip().split(' ')
        value = [float(i) for i in value]
        x.append(value[0])
        y.append(value[1])
        z.append(value[2])


n = 5

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
ax.set_xlim(0, n)
ax.set_ylim(0, n)
ax.set_zlim(0, n)

plt.show()
