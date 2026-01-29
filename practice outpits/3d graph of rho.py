from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import numpy as np

ax = plt.figure().add_subplot(projection = '3d')

def my_funct(x,y):
    return x*y

x = np.linspace(0, 100, 100)
y = np.linspace(0,1,100)
X, Y = np.meshgrid(x, y)

zs = np.array(my_funct(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.legend()
ax.set_xlim(0,100)
ax.set_ylim(0,1)
ax.set_zlim(0,100)
ax.set_xlabel('Number in queue')
ax.set_ylabel('Individual reneging rate')

plt.show()