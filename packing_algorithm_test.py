import numpy as np
import matplotlib.pyplot as plt
import math
#from mpl_toolkits.mplot3d import Axes3D

ball=6
data = np.random.rand(3,ball)*255

x, y, z = data[0], data[1], data[2]
print(x,y,z)
ax = plt.subplot(111, projection='3d')  # 创建一个三维的绘图工程
#  将数据点分成三部分画，在颜色上有区分度
ax.scatter(x, y, z, c='r')  # 绘制数据点

ax.set_zlabel('Z')  # 坐标轴
ax.set_ylabel('Y')
ax.set_xlabel('X')
plt.show()

learningRate=0.01

for ii in range(1,1000000):
  for jj in range(0,ball):
    sumx=sum(x)
    sumy=sum(y)
    sumz=sum(z)
    GradX=ball*x[jj]-sumx
    GradY=ball*y[jj]-sumy
    GradZ=ball*z[jj]-sumz
    tempsum=GradX*GradX+GradY*GradY+GradZ*GradZ
    tempsum=math.sqrt(tempsum)


    if (ii==10 or ii==100 or ii==1000 or ii==10000 or ii==100000) and jj==0:
      print('tempsum:',tempsum,'GradX:',GradX)

    GradX2=GradX/tempsum
    GradY2=GradY/tempsum
    GradZ2=GradZ/tempsum
    x[jj]=x[jj]-GradX2
    y[jj]=y[jj]-GradY2
    z[jj]=z[jj]-GradZ2

  
  if ii==10 or ii==100 or ii==1000 or ii==10000 or ii==100000:



    #print(x,y,z)
    ax = plt.subplot(111, projection='3d')  # 创建一个三维的绘图工程
    #  将数据点分成三部分画，在颜色上有区分度
    ax.scatter(x, y, z, c='r')  # 绘制数据点

    ax.set_zlabel('Z')  # 坐标轴
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()


    
