# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 14:46:05 2019

@author: ws
"""
from cast_heat import two_densonal_diff
from cast_heat import steady_temp_cal
from crack import crack
from stress import solver

var_h_initial=[245.4,181.4,146.4,182.5] # 待定参数的初值# 待定参数的初值
def ObjectiveFunction(var_h_initial):
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
    #var_aimTemp=[800,750,700,650] # 测量的温度值 # 测量的温度值
    #var_SCZ_num=4 # 连铸二冷区冷却区段的数量
    #Var_runningSteadyTime = 5 #稳态运行时间
    ########## 需要输入的参数 end   #################
    #
    ##############  parameters for optimation start ################
    var_deltTime=0.4 # 差分计算时间间隔
    #var_iter_max=3; # 优化迭代的最大次数
    #var_eps=200 # 优化的停止标准
    
    #############  parameters for optimation end ################
    
    MiddleTemp=[([var_castingTemp]*var_YNumber) for i in range(var_XNumber)]
    tl=[0]*len(var_dis);  # 铸坯凝固时间的初值var_dis=[0.0,0.9,1.27,3.12,5.32,7.68] # 连铸二冷区各段距弯月面的距离 单位m
    for i in range(len(var_dis)):
        tl[i]=var_dis[i]/var_VcastOriginal#var_VcastOriginal=0.6/60
    time_Mold = int((tl[1]-tl[0])/var_deltTime)
    time_SCZ = int((tl[len(var_dis)-1]-tl[1])/var_deltTime)#var_deltTime=0.4 # 差分计算时间间隔
    Time_all = time_Mold+time_SCZ
    MiddleTemp_all22,t= steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial,var_ZNumber,var_castingTemp)
    result_crack=crack(MiddleTemp_all22,var_XNumber,var_X,var_YNumber,var_Y,Time_all,var_liqTemp,var_SodTemp)#缩孔
    stress,estress,strain=solver(20,20,0.1,0.1,t,15,1900)#计算节点位移
    return strain,result_crack
    
print(ObjectiveFunction(var_h_initial))

a,b=ObjectiveFunction(var_h_initial)
print(a[-2][-2])



















