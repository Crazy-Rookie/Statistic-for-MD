#################################################################################
#To compute standard deviation of POSCAR (p0001, p0002,...)
#e.g.: python std_coor.py 1 2000
#可以计算原子运动速度与位置标准差
#2021.12.21完成，终极版
#202021150066@mail.bnu.edu.cn
#################################################################################

import sys
import math
import copy
import numpy as np

first_no=int(sys.argv[1])
last_no=int(sys.argv[2])

file1=open('p%04d' %first_no, 'r')
lines1=file1.readlines()
file1.close()
no_atom=int(0)
no_atom1=lines1[6].split()
for i in range(len(no_atom1)):
    no_atom=no_atom+int(no_atom1[i])

lattice = np.zeros((3, 3))
for i in range(2, 5):
    for j in range(3):
        lattice[i - 2][j] = float(lines1[i].split()[j])    
    
distance=[0]*no_atom
std_coor=[0]*no_atom
ver=[0]*no_atom
std=open('std.txt', 'w')
coor=open('avg_coor.txt', 'w')
verlocity=open('avg_verlocity.txt', 'w')


coor0=[]
coor1=[]
coor2=[]
dir_coor0=[[0 for col in range(3)] for row in range(no_atom)]
dir_coor1=[[0 for col in range(3)] for row in range(no_atom)]
dir_coor2=[[0 for col in range(3)] for row in range(no_atom)]
dis_coor=[[0 for col in range(3)] for row in range(no_atom)]
avg_coor=np.array([[0 for col in range(3)] for row in range(no_atom)])
avg_ver=np.array([[0 for col in range(3)] for row in range(no_atom)])



keywords=['D', 'd','S','s']
#print lines1[6]
#print lines1[6][0]
file0=open('POSCAR', 'r')
lines0=file0.readlines()
file0.close()

for m in range (0, no_atom):
    coor0=lines0[m+8].split()
    for z in range(3):
        dir_coor0[m][z]=float(coor0[z])
        
for i in  range (first_no, last_no):
    file1=open('p%04d' %i, 'r')
    lines1=file1.readlines()
    file1.close()
    
    file2=open('p%04d' %(i+1), 'r')
    lines2=file2.readlines()
    file2.close()
   
    for m in range (0, no_atom):
        coor1=lines1[m+8].split()
        coor2=lines2[m+8].split()
        for z in range(3):
            dir_coor1[m][z]=float(coor1[z])
            dir_coor2[m][z]=float(coor2[z])
        
        
            if (dir_coor1[m][z]-dir_coor0[m][z])>0.5:
                dir_coor1[m][z]=dir_coor1[m][z]-1
            if (dir_coor1[m][z]-dir_coor0[m][z])<-0.5:
                dir_coor1[m][z]=dir_coor1[m][z]+1    
            if (dir_coor2[m][z]-dir_coor0[m][z])>0.5: 
                dir_coor2[m][z]=dir_coor2[m][z]-1    
            if (dir_coor2[m][z]-dir_coor0[m][z])<-0.5:    
                dir_coor2[m][z]=dir_coor2[m][z]+1        
        
    delta_coor=np.array(dir_coor2)-np.array(dir_coor1)
    delta_coor=np.maximum(delta_coor,-delta_coor)                  #verlocity need to be positive
    avg_ver=avg_ver+delta_coor
    avg_coor=avg_coor+np.array(dir_coor1)    
        
avg_ver=avg_ver/(last_no-first_no)
for m in range(no_atom):
    ver[m]=(np.sum(np.dot(avg_ver[m],lattice)**2))**0.5             #verlocity saved file
avg_coor=avg_coor/(last_no-first_no)                                #coor saved file



dir_coor1=[[0 for col in range(3)] for row in range(no_atom)]

for i in  range (first_no+1, last_no):
    file1=open('p%04d' %i, 'r')
    lines1=file1.readlines()
    file1.close()  
    for m in range (0, no_atom):
        coor1=lines1[m+8].split()
        for z in range(3):
            dir_coor1[m][z]=float(coor1[z])
        

            if (dir_coor1[m][z]-dir_coor0[m][z])>0.5:
                dir_coor1[m][z]=dir_coor1[m][z]-1
            if (dir_coor1[m][z]-dir_coor0[m][z])<-0.5:
                dir_coor1[m][z]=dir_coor1[m][z]+1                   
        #this loop only for std

        delta2_coor=np.array(dir_coor1[m])-avg_coor[m]
        distance[m] =distance[m]+ np.sum(np.dot(delta2_coor,lattice)**2)  

    

dis=(np.array(distance)/(last_no-first_no))**0.5           #std saved file

for m in range(no_atom):
    for z in range(3):
        print(avg_coor[m][z], end=' ', file=coor)
    print(dis[m],file=std)
    print(file=coor)
    print(ver[m],file=verlocity)





















