# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 16:40:10 2019

@author: ws
"""
from cast_heat import two_densonal_diff
from cast_heat import steady_temp_cal
from crack import crack
import numpy as np
import math
var_ZNumber = 500 # 拉坯方向网格划分的数量
var_XNumber=20 # 铸坯厚度方向的网格点数
var_X=0.16 # 铸坯厚度
var_YNumber=20 # 铸坯厚度方向的网格点数
var_Y=0.16 # 铸坯厚度
var_Z=7.68
deltZ = var_Z/var_ZNumber # 拉皮方向的空间间隔
middle_temp = [([1530]*var_XNumber) for i in range(var_XNumber)]
var_temperatureWater = 30
var_rouS=7800 
var_rouL=7200 
var_specificHeatS=660 
var_specificHeatL=830 
var_TconductivityS=31 
var_TconductivityL=35 
var_liqTemp=1514 
var_SodTemp=1453 
var_m=7.1
var_latentHeatofSolidification = 268000
var_VcastOriginal=0.6/60
var_controlTime=4
var_dis=[0.0,0.9,1.27,3.12,5.32,7.68] # 连铸二冷区各段距弯月面的距离 单位m
var_runningTime=20
var_castingTemp=1530 # 铸坯的浇筑温度
var_h_initial=[245.4,181.4,146.4,182.5] # 待定参数的初值# 待定参数的初值
#############  parameters for optimation end ################
var_deltTime=0.4 # 差分计算时间间隔
MiddleTemp=[([var_castingTemp]*var_YNumber) for i in range(var_XNumber)]

tl=[0]*len(var_dis);  # 铸坯凝固时间的初值var_dis=[0.0,0.9,1.27,3.12,5.32,7.68] # 连铸二冷区各段距弯月面的距离 单位m
for i in range(len(var_dis)):
    tl[i]=var_dis[i]/var_VcastOriginal#var_VcastOriginal=0.6/60
time_Mold = int((tl[1]-tl[0])/var_deltTime)
time_SCZ = int((tl[len(var_dis)-1]-tl[1])/var_deltTime)#var_deltTime=0.4 # 差分计算时间间隔
Time_all = time_Mold+time_SCZ
steady_temp_cal,t= steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial,var_ZNumber,var_castingTemp)


def lnd_xy_creack(nx,ny,dx,dy):#创建节点编码矩阵nx,ny为边界的节点个数dx为单元边长
    unit=(nx-1)*(ny-1)*2#单元个数
    lnd= [[0 for i in range(3)] for j in range(unit)]
    l=0
    for j in range(ny-1):
        for i in range(nx-1):
            lnd[l][0]=i+j*(nx)
            lnd[l][1]=i+1+j*(nx)
            lnd[l][2]=i+(j+1)*(nx)
            lnd[l+1][0]=i+j*(nx)+1
            lnd[l+1][1]=i+1+(j+1)*(nx)
            lnd[l+1][2]=i+(j+1)*(nx)
            l=l+2
    m=nx*ny#总节点数
    xy=[[0 for i in range(4)] for j in range(m)]
     
    for i in range(ny):
        l=i*nx
        for j in range(nx):
            xy[l+j][0]=j*dx
            xy[l+j][1]=i*dy
            xy[l+j][2]=j#x轴坐标
            xy[l+j][3]=i#y轴坐标
    return lnd,xy,unit,m
#几个计算过程用到的变量
def Ee(t):
    if t>900:
        e = 31.46 - 0.02256 * t + 0.00000138 * pow(t, 2)
    else:
        e = 347.6526 - 0.3727 * t
    #e=e * 1000000000.0
    return e
def u(t):
    u = 0.287 + 0.0000823 * t
    return u
def a(t):
    aa = (0.435 + 0.000708 * t) * 0.00001
    return aa
def area(ie,lnd,xy):
    i=lnd[ie][0]#i.j.k代表节点号
    j=lnd[ie][1]
    k=lnd[ie][2]
    area = (xy[j][0] - xy[i][0]) * (xy[k][1] - xy[i][1] ) - (xy[k][0]  - xy[i][0]) * (xy[j][1] - xy[i][1])
    area = area * 0.5
    return area
def betaf(ie,t):#这部分有疑问待解决
    beta=(1-2*u(t))/(1+u(t))
    return beta
    
#计算弹朔性矩阵
def md(ie,t):#这里的t是温度单个单元的温度
    beta = betaf(ie, t)
    d=[[0 for i in range(3)] for j in range(3)]
    d[0][0]=2+beta
    d[0][1]=1-beta
    d[1][0]=d[0][1]
    d[1][1]=d[0][0]
    d[2][2]=0.5+beta
    for i in range(3):
        for j in range(3):
            d[i][j]=beta*Ee(t)/(1-2*u(t))/(1+2*beta)*d[i][j]

    return d        

#计算位移——应变转换矩阵（几何矩阵B）
def mb(ie,lnd,xy):#ie代表的是单元序号
    i=lnd[ie][0]#i.j.k代表节点号
    j=lnd[ie][1]
    k=lnd[ie][2]
    b=[[0 for i in range(6)] for j in range(3)]
    b[0][0]=xy[j][1]-xy[k][1]
    b[0][2]=xy[k][1]-xy[i][1]
    b[0][4]=xy[i][1]-xy[j][1]
    b[1][1]=xy[k][0]-xy[j][0]
    b[1][3]=xy[i][0]-xy[k][0]
    b[1][5]=xy[j][0]-xy[i][0]
    b[2][0]= b[1][1]
    b[2][1]= b[0][0]
    b[2][2]= b[1][3]
    b[2][3]= b[0][2]
    b[2][4]= b[1][5]
    b[2][5]= b[0][4]
    s=area(ie,lnd,xy)
    #s=0.5
    for ii in range(3):
        for jj in range(6):
            b[ii][jj]=b[ii][jj]/(2*s)
    return b
#求解单元刚度矩阵
def mke(ie,t,lnd,xy,kk):#kk代表当前时间切片,这里的t是温度场
    i=lnd[ie][0]#i.j.k代表节点号
    j=lnd[ie][1]
    k=lnd[ie][2]
    tpkk=t[kk]
    te=tpkk[xy[i][2]][xy[i][3]]
    te+=tpkk[xy[j][2]][xy[j][3]]
    te+=tpkk[xy[k][2]][xy[k][3]]
    te=te/3.0
    d=md(ie,te)
    b=mb(ie,lnd,xy)#3*6的矩阵
    s=[[0 for i in range(6)] for j in range(3)]
    
    for ii in range(3):
        for jj in range(6):
            s[ii][jj]=0
            for kk in range(3):
                s[ii][jj]=s[ii][jj]+d[ii][kk]*b[kk][ii]
    ake=[[0 for i in range(6)] for j in range(6)]
    for ii in range(6):
        for jj in range(6):
            ake[ii][jj]=0
            for kk in range(3):
                 ake[ii][jj]= ake[ii][jj]+s[kk][ii]*b[kk][jj]            
    return ake       
#求解整体刚度矩阵
def mk(unit,m,t,lnd,xy,kk):
    N=2*m
    ak=[[0 for i in range(N)] for j in range(N)]
    for ie in range(unit):
        ake=mke(ie,t,lnd,xy,kk)
        i=lnd[ie][0]#i.j.k代表节点号
        j=lnd[ie][1]
        k=lnd[ie][2]
        #print(i,j,k)
        ak[2*i][2*i]=ak[2*i][2*i]+ake[0][0]
        ak[2*i][2*i+1]+=ake[0][1]
        ak[2*i][2*j]+=ake[0][2]
        ak[2*i][2*j+1]+=ake[0][3]
        ak[2*i][2*k]+=ake[0][4]
        ak[2*i][2*k+1]+=ake[0][5]
        ak[2*i+1][2*i]+=ake[1][0]
        ak[2*i+1][2*i+1]+=ake[1][1]
        ak[2*i+1][2*j]+=ake[1][2]
        ak[2*i+1][2*j+1]+=ake[1][3]
        ak[2*i+1][2*k]+=ake[1][4]
        ak[2*i+1][2*k+1]+=ake[1][5]
        
        ak[2*j][2*i]+=ake[2][0]
        ak[2*j][2*i+1]+=ake[2][1]
        ak[2*j][2*j]+=ake[2][2]
        ak[2*j][2*j+1]+=ake[2][3]
        ak[2*j][2*k]+=ake[2][4]
        ak[2*j][2*k+1]+=ake[2][5]
        ak[2*j+1][2*i]+=ake[3][0]
        ak[2*j+1][2*i+1]+=ake[3][1]
        ak[2*j+1][2*j]+=ake[3][2]
        ak[2*j+1][2*j+1]+=ake[3][3]
        ak[2*j+1][2*k]+=ake[3][4]
        ak[2*j+1][2*k+1]+=ake[3][5]
        
        ak[2*k][2*i]+=ake[4][0]
        ak[2*k][2*i+1]+=ake[4][1]
        ak[2*k][2*j]+=ake[4][2]
        ak[2*k][2*j+1]+=ake[4][3]
        ak[2*k][2*k]+=ake[4][4]
        ak[2*k][2*k+1]+=ake[4][5]
        ak[2*k+1][2*i]+=ake[5][0]
        ak[2*k+1][2*i+1]+=ake[5][1]
        ak[2*k+1][2*j]+=ake[5][2]
        ak[2*k+1][2*j+1]+=ake[5][3]
        ak[2*k+1][2*k]+=ake[5][4]
        ak[2*k+1][2*k+1]+=ake[5][5]
    return ak

#计算热荷载
def heatload(t,kk,lnd,xy,unit,kk1,m):#t,kk为温度场和时间切片,kk1为段末坐标
    tpkk=t[kk]
    tpkk1=t[kk1]
    p=[0 for i in range(2*m)]
    for ie in range (unit):
        i=lnd[ie][0]#i.j.k代表节点号
        j=lnd[ie][1]
        k=lnd[ie][2]
        te=tpkk[xy[i][2]][xy[i][3]]
        te+=tpkk[xy[j][2]][xy[j][3]]
        te+=tpkk[xy[k][2]][xy[k][3]]
        te=te/3.0
        te1=tpkk1[xy[i][2]][xy[i][3]]
        te1+=tpkk1[xy[j][2]][xy[j][3]]
        te1+=tpkk1[xy[k][2]][xy[k][3]]
        te1=te/3.0
        beta = betaf(ie, te)
        e=Ee(te)
        aa=a(te)
        uu=u(te)
        load=3*beta*e*aa*(te-te1)/(2*(1-2*uu*(1+2*beta)))
        p[2*i]=p[2*i]+load*(xy[j][1]-xy[k][1])
        p[2*i+1]=p[2*i+1]+load*(xy[k][0]-xy[j][0])
        p[2*j]=p[2*j]+load*(xy[k][1]-xy[i][1])
        p[2*j+1]=p[2*j+1]+load*(xy[i][0]-xy[k][0])
        p[2*k]=p[2*k]+load*(xy[i][1]-xy[j][1])
        p[2*k+1]=p[2*k+1]+load*(xy[j][0]-xy[i][0])
    
    return p
def strain_stress(t,strain):
    if t<1100:#貌似没啥子用
        estrain = 0.000484 - 0.000000368 * t
    else:
        estrain = 0.000147 - 0.00000008 * t
    estrain = 0.00161 - 0.000000756 * t#这个有用
    if strain < estrain:
        stress = Ee(t) * strain
    else:
        stress = Ee(t) * estrain + (strain - estrain) * 0.1 * Ee(t)
    return stress          
    
#求解节点位移
def solver(nx,ny,dx,dy,t,kk,kk1):
    lnd,xy,unit,m=lnd_xy_creack(nx,ny,dx,dy)
    mk1=mk(unit,m,t,lnd,xy,kk)
    pp=heatload(t,kk,lnd,xy,unit,kk1,m)
    k=np.array(mk1)
    #print(k)
    p=np.array(pp)
    #求解节点位移
    wx=(np.linalg.pinv(k)).dot(p.reshape(p.shape[0],1))
    tpkk=t[kk]
    tpkk1=t[kk1]
    xe=[0 for i in range(6)]#创建节点位移
    stress=[0 for i in range(unit)]#单元应变
    strain=[0 for i in range(unit)]#单元应力
    estress=[[0 for i in range(2)] for j in range(unit)]
    print(unit)
    for ie in range (unit):
        i=lnd[ie][0]#i.j.k代表节点号
        j=lnd[ie][1]
        k=lnd[ie][2]
        te=tpkk[xy[i][2]][xy[i][3]]
        te+=tpkk[xy[j][2]][xy[j][3]]
        te+=tpkk[xy[k][2]][xy[k][3]]
        te=te/3.0
        te1=tpkk1[xy[i][2]][xy[i][3]]
        te1+=tpkk1[xy[j][2]][xy[j][3]]
        te1+=tpkk1[xy[k][2]][xy[k][3]]
        te1=te/3.0
        b=mb(ie,lnd,xy)#3*6的矩阵
        xe[0]=wx[2*i]
        xe[1]=wx[2*i+1]
        xe[2]=wx[2*j]
        xe[3]=wx[2*j+1]
        xe[4]=wx[2*k]
        xe[5]=wx[2*k+1]
        xe=np.array(xe)
        stress[ie]=np.array(b).dot(xe.reshape(xe.shape[0],1))
        #beta = betaf(ie, te)
        ybx=stress[ie][0]#x方向的应变
        yby=stress[ie][1]#y方向的应变
        ybxy=stress[ie][2]#切方向的应变
        #ybx=-1*(1-beta)/(1+2*beta)*(ybx+yby)#z方向的应变
        ybz=-u(te)/(1-u(te))*(ybx+yby)#z方向的应变
        #等效应变
        estress[ie][0]=math.sqrt(math.pow(ybx-yby,2)+ math.pow(yby-ybz,2) + math.pow(ybz-ybx,2) + 1.5 * math.pow(ybxy,2))/ (1 + u(te)) / math.sqrt(2)
        #等效应力
        estress[ie][1] = strain_stress(te,estress[ie][0])
        d=md(ie,te)
        strain[ie]=np.array(d).dot(stress[ie])
    return stress,estress,strain
    
stress,estress,strain=solver(20,20,0.1,0.1,t,1500,1900)#计算节点位移
print(stress[500],estress[500],strain[500])



















