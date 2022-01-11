'''
这个脚本可以用于统计MD输出结构的键长和键角信息。
脚本中需要特殊设定的数值均列在开头处，并用注释写明了每个数值设置的含义
注意：该脚本仅适用于钙钛矿结构，其他结构可能也能给出结果，但统计效果不佳。
lhr：202021150066@mail.bnu.edu.cn

输出文件：
bondst.dat: A-B bond 的键长分布情况
anglest.dat: A-B-A angle 的键角分布情况，其中B为top atom
dihedralst.dat：B-A--A-B 即两个近邻A原子指向同一方向的A-B键之间的二面角

'''

import sys
import math
import copy
import numpy as np

first_no = 1                           # which p-file do you want to start
last_no = 8000                        # which to finish

atoms1 = list(range(0, 8))             # index of A atom(s) in A-B bond
atoms2 = list(range(16, 40))          # index of B atom(s) in A-B bond
Maxbond = 4.2                           # This value should be a little biger than maximum length of A-B bonds.
Minbond = 2

print('The script is running, it will take a few minutes.')
Distance = np.zeros((last_no-first_no+1,1))
bondst=[]
anglest=[]
dihedralst=[]
dis=[]


flb = open("bond_dis.dat", "w")
fla = open("ang_ang.dat", "w")

file1 = open('POSCAR', 'r')
lines1 = file1.readlines()
file1.close()
lattice = np.zeros((3, 3))
for i in range(2, 5):
    for j in range(3):
        lattice[i - 2][j] = float(lines1[i].split()[j])



no_atom = sum(int(i) for i in lines1[6].split())
distance = [0] * no_atom
avg_velocity = [0] * no_atom

keywords = ['D', 'd', 'C', 'c']

file2 = open('POSCAR' , 'r')
lines2 = file2.readlines()
file2.close()
dis=[]
cart_coor1 = [[0 for col in range(3)] for row in range(no_atom)]
bonddata=[[] for row in range(no_atom)]


if lines2[7][0] in keywords:
    for m in range(no_atom):
        direct_coor = lines2[m + 8].split()
        for z in range(3):
            cart_coor1[m][z] = float(direct_coor[z])
else:
    print('Plase prepare Direct coordinate.')

for j in atoms1:
    for k in atoms2:
        length = np.array(cart_coor1[j]) - np.array(cart_coor1[k])
        for l in range(3):
            if length[l]>0.5:
                length[l]=length[l]-1
            if length[l]<-0.5:
                length[l]=length[l]+1    
           
        length = (np.sum(np.dot(length,lattice)**2))**0.5
        
        if length < Maxbond and length > Minbond:           
            dis.append(length)
            bonddata[j].append(k)
 
for i in atoms1:
    for j in atoms1  : 
        if j > i: 
            connectatom=[x for x in bonddata[i] if x in bonddata[j]] 
            if len(connectatom) !=0 : 
                 
                vectora=np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[connectatom[0]])),lattice) 
                for l in range(3):
                    if vectora[l]>0.5:
                        vectora[l]=vectora[l]-1
                    if vectora[l]<-0.5:
                        vectora[l]=vectora[l]+1                    
                vectora=vectora/np.linalg.norm(vectora) 
                vectorb=np.dot((np.array(cart_coor1[j])-np.array(cart_coor1[connectatom[0]])),lattice) 
                for l in range(3):
                    if vectorb[l]>0.5:
                        vectorb[l]=vectorb[l]-1
                    if vectorb[l]<-0.5:
                        vectorb[l]=vectorb[l]+1                                      
                vectorb=vectorb/np.linalg.norm(vectorb) 
                cos_angle=np.sum(vectora*vectorb)  #这个值是python运算时产生的误差，无法消除            
                angle=np.arccos(cos_angle)*180/math.pi            
                anglest.append(angle)   
                vectora=np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[connectatom[1]])),lattice) 
                for l in range(3):
                    if vectora[l]>0.5:
                        vectora[l]=vectora[l]-1
                    if vectora[l]<-0.5:
                        vectora[l]=vectora[l]+1      
                vectora=vectora/np.linalg.norm(vectora) 
                vectorb=np.dot((np.array(cart_coor1[j])-np.array(cart_coor1[connectatom[1]])),lattice) 
                for l in range(3):
                    if vectorb[l]>0.5:
                        vectorb[l]=vectorb[l]-1
                    if vectorb[l]<-0.5:
                        vectorb[l]=vectorb[l]+1                   
                vectorb=vectorb/np.linalg.norm(vectorb) 
                cos_angle=np.sum(vectora*vectorb)  #这个值是python运算时产生的误差，无法消除            
                angle=np.arccos(cos_angle)*180/math.pi            
                anglest.append(angle)                   
                
                               
no_bond=len(dis)
no_ang=len(anglest)
bond_evo=[[] for row in range(last_no-first_no+1)]
ang_evo=[[] for row in range(last_no-first_no+1)]




for p in range(first_no , last_no+1):

    file2 = open('p%04d' % p, 'r')
    lines2 = file2.readlines()
    file2.close()
    anglest=[]
    cart_coor1 = [[0 for col in range(3)] for row in range(no_atom)]

    
    if lines2[7][0] in keywords:
        for m in range(no_atom):
            direct_coor = lines2[m + 8].split()
            for z in range(3):
                cart_coor1[m][z] = float(direct_coor[z])
    else:
        print('Plase prepare Direct coordinate.')

    for j in atoms1:
        for k in bonddata[j]:
            length = np.array(cart_coor1[j]) - np.array(cart_coor1[k])
            for l in range(3):
                if length[l]>0.5:
                    length[l]=length[l]-1
                if length[l]<-0.5:
                    length[l]=length[l]+1               
            length = (np.sum(np.dot(length,lattice)**2))**0.5
            
                    
            dis.append(length)
            bond_evo[p-1].append(length)
    bondst.extend(dis)

    
    # The statistics of distance is done above. 
    
    
    
    
    
    for i in atoms1:
        for j in atoms1  : 
            if j > i: 
                connectatom=[x for x in bonddata[i] if x in bonddata[j]] 
                if len(connectatom) !=0 : 
                     
                    vectora=np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[connectatom[0]])),lattice) 
                    for l in range(3):
                        if vectora[l]>0.5:
                            vectora[l]=vectora[l]-1
                        if vectora[l]<-0.5:
                            vectora[l]=vectora[l]+1                    
                    vectora=vectora/np.linalg.norm(vectora) 
                    vectorb=np.dot((np.array(cart_coor1[j])-np.array(cart_coor1[connectatom[0]])),lattice) 
                    for l in range(3):
                        if vectorb[l]>0.5:
                            vectorb[l]=vectorb[l]-1                
                        if vectorb[l]<-0.5:                
                            vectorb[l]=vectorb[l]+1                                                      
                    vectorb=vectorb/np.linalg.norm(vectorb)                 
                    cos_angle=np.sum(vectora*vectorb)  #这个值是python运算时产生的误差，无法消除                            
                    angle=np.arccos(cos_angle)*180/math.pi   
                    if angle<90:
                        angle=180-angle                    
                    anglest.append(angle)     
                    
                    ang_evo[p-1].append(angle)
                    
                    vectora=np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[connectatom[1]])),lattice)                 
                    for l in range(3):                
                        if vectora[l]>0.5:                
                            vectora[l]=vectora[l]-1                
                        if vectora[l]<-0.5:                
                            vectora[l]=vectora[l]+1                      
                    vectora=vectora/np.linalg.norm(vectora)                 
                    vectorb=np.dot((np.array(cart_coor1[j])-np.array(cart_coor1[connectatom[1]])),lattice)                 
                    for l in range(3):                
                        if vectorb[l]>0.5:                
                            vectorb[l]=vectorb[l]-1                
                        if vectorb[l]<-0.5:                
                            vectorb[l]=vectorb[l]+1                                   
                    vectorb=vectorb/np.linalg.norm(vectorb)                 
                    cos_angle=np.sum(vectora*vectorb)  #这个值是python运算时产生的误差，无法消除                            
                    angle=np.arccos(cos_angle)*180/math.pi                            
                    if angle<90:
                        angle=180-angle
                    anglest.append(angle)                                   
                                                      
                    ang_evo[p-1].append(angle)
                    


                
    # The statistics of angle is done above. 

#
#         
#            
#                    for k in bonddata[i] :
#                        if k != connectatom[0]:
#                            vector1= np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[k])),lattice)
#                            vector1=vector1/np.linalg.norm(vector1)
#                            for l in bonddata[j] :
#                                if l!= connectatom[0]:
#                                    vector2= np.dot((np.array(cart_coor1[j])-np.array(cart_coor1[l])),lattice)
#                                    vector2=vector2/np.linalg.norm(vector2)
#                                    vector3= np.dot((np.array(cart_coor1[i])-np.array(cart_coor1[j])),lattice)
#                                    vector3=vector3/np.linalg.norm(vector3)
#                                    cos_angle=np.sum(vector1*vector2)
#                                    angle=np.arccos(cos_angle)*180/math.pi
#                                    
#                                    if angle <45 and angle>-45:
#                                        vector4=np.cross(vector1,vector3)
#                                        vector5=np.cross(vector2,vector3)
#                                        vector4=vector4/np.linalg.norm(vector4)
#                                        vector5=vector5/np.linalg.norm(vector5)
#                                        cos_dihedral=np.sum(vector4*vector5)
#                                        dihedral=np.arccos(cos_dihedral)*180/math.pi
#                                        dihedralst.append(dihedral)
                                    
    # The statistics of dihedral is done above. 
                    
                
    
for j in range(last_no-first_no+1):
    for k in range(no_bond):
        print(bond_evo[j][k], end='       ', file=flb)
    print(file=flb)   
    
for j in range(last_no-first_no+1):
    for k in range(no_ang):
        print(ang_evo[j][k], end='       ', file=fla)
    print(file=fla)       


#datafl = open('bondst.dat','w')
#for i in range(int(Maxbond*100+1)):
#    count=len([x for x in bondst if x>=(i-1)/100 and x<i/100])
#    print(str(i/100)+'           '+str(count/len(bondst)),file=datafl)
    
#datafl = open('anglest.dat','w')
#for i in range(180*10+1):
#    if len(anglest)==0:
#        print('No angle available.')
#    else:
#        count=len([x for x in anglest if x>(i-1)/10 and x<=i/10])
#        print(str(i/10)+'           '+str(count/len(anglest)),file=datafl)

#datafl = open('dihedralst.dat','w')
#if len(dihedralst) == 0:
#    print('No diheral available.')
#else:
#    for i in range(180*10+1):
#            count=len([x for x in dihedralst if x>=(i-1)/10 and x<i/10])
#            print(str(i/10)+'           '+str(count/len(dihedralst)),file=datafl)





 
    


    
