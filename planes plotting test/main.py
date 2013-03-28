import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def fun(x, y):
  return 20.0 - x - y

def fun2(x,y):
  return 22.0 - x

def fun3(x,y):
  return 25.0 - y

def fun4(x,y):
  return 18.0 - y

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 20
xs = [i for i in range(n) for _ in range(n)]
ys = range(n) * n
zs = [fun(x, y) for x,y in zip(xs,ys)]
ax.scatter(xs, ys, zs, c='r')

xs = [i for i in range(n) for _ in range(n)]
ys = range(n) * n
zs = [fun2(x, y) for x,y in zip(xs,ys)]
ax.scatter(xs, ys, zs)

xs = [i for i in range(n) for _ in range(n)]
ys = range(n) * n
zs = [fun3(x, y) for x,y in zip(xs,ys)]
ax.scatter(xs, ys, zs,c ='g')

xs = [i for i in range(n) for _ in range(n)]
ys = range(n) * n
zs = [fun4(x, y) for x,y in zip(xs,ys)]
ax.scatter(xs, ys, zs,c ='y')


ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
