import pdb2ball_single
import numpy as np
import matplotlib.pyplot as plt
import math
import Bio



def get_random_radius_list():
    '''get radius value randomly'''
    #p_d = pdb2ball_single.pdb2ball_sinble(".\\PDB_ori\\")
    p_d = pdb2ball_single.pdb2ball_sinble(".\\pdbfile\\")
    r_s = []
    for protein in p_d:
        r_s.append(p_d[protein]['radius'])
    rand_r_s = []
    rand_list = np.random.randint(len(r_s), size=ball)
    for i in rand_list:
        rand_r_s.append(r_s[i])
    return rand_r_s

#ball_radius = get_random_radius_list()

#print('radius:',ball_radius)
################以上是获得半径

def simple_packing_func():
    '''pack 5 balls with a simple algorithm'''
    ball_radius = get_random_radius_list()
    print('radius:',ball_radius)

#设置球的数量为5
    ball=5

####一直随机扔球直到初始没有重叠的球
    while 1:
        mark = 0
        data = np.random.rand(3,ball)*25500
        x, y, z = data[0], data[1], data[2]
        for ii in range(0,ball):
            for jj in range(ii+1,ball):
                tempR1 = ball_radius[ii]
                tempR2 = ball_radius[jj]
                tempDis = (x[ii]-x[jj])**2 + (y[ii]-y[jj])**2 + (z[ii]-z[jj])**2
                Distance = math.sqrt(tempDis)
                if tempR1 + tempR2 > Distance:
                    mark = 1
                    print('放置的球重叠了！球的半径以及球之间的距离分别是：',tempR1,tempR2,Distance)

        if mark == 0:
            break


    print('Finish rolling balls')

    print('Coordinates:')
    print(x,y,z)
    ax = plt.subplot(111, projection='3d')  # 创建一个三维的绘图工程
#  将数据点分成三部分画，在颜色上有区分度
    ax.scatter(x, y, z, c='r')  # 绘制数据点

    ax.set_zlabel('Z')  # 坐标轴
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()

#learningRate=0.01

    for ii in range(1,100001):
        for jj in range(0,ball):
            sumx=sum(x)
            sumy=sum(y)
            sumz=sum(z)
            GradX=ball*x[jj]-sumx
            GradY=ball*y[jj]-sumy
            GradZ=ball*z[jj]-sumz
            tempsum=GradX*GradX+GradY*GradY+GradZ*GradZ
            tempsum=math.sqrt(tempsum) 


            if (ii==10 or ii==100 or ii==1000 or ii==10000 or ii==100000 or ii==1000000) and jj==0:
                print('tempsum:',tempsum,'GradX:',GradX)

            GradX2=GradX/tempsum
            GradY2=GradY/tempsum
            GradZ2=GradZ/tempsum
            x[jj]=x[jj]-GradX2
            y[jj]=y[jj]-GradY2
            z[jj]=z[jj]-GradZ2
######原代码基础上加入判断每次移动后是否会导致重叠
            mark = 0
            for kk in range(0,ball):
                for ll in range(kk+1,ball):
                    tempR1 = ball_radius[kk]
                    tempR2 = ball_radius[ll]
                    tempDis = (x[kk]-x[ll])**2 + (y[kk]-y[ll])**2 + (z[kk]-z[ll])**2
                    Distance = math.sqrt(tempDis)
                    Distance = math.sqrt(tempDis)
                    if tempR1 + tempR2 > Distance:
                        mark = 1

            if mark == 1:
                x[jj]=x[jj]+GradX2
                y[jj]=y[jj]+GradY2
                z[jj]=z[jj]+GradZ2
############

############
  
        if ii==10 or ii==100 or ii==1000 or ii==10000 or ii==100000 or ii==1000000:
            ax = plt.subplot(111, projection='3d')  # 创建一个三维的绘图工程
            #将数据点分成三部分画，在颜色上有区分度
            ax.scatter(x, y, z, c='r')  # 绘制数据点

            ax.set_zlabel('Z')  # 坐标轴
            ax.set_ylabel('Y')
            ax.set_xlabel('X')
            plt.show()



########每两个点之间的距离
    for ii in range(0,ball):
        for jj in range(ii+1,ball):
            tempR1 = ball_radius[ii]
            tempR2 = ball_radius[jj]
            tempDis = (x[ii]-x[jj])**2 + (y[ii]-y[jj])**2 + (z[ii]-z[jj])**2
            Distance = math.sqrt(tempDis)
            print('半径：',ii,'球',tempR1,'半径：',jj,'球',tempR2,'距离：',Distance)
            
            



        

 
