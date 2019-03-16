from two_densonal_diff_lifecycle_for_mpc import two_densonal_lifecycle_for_mpc
from two_densonal_diff_lifecycle_for_mpc import steady_temp_cal
import numpy.linalg as nlg
import numpy as np
def initialize_parameters_for_continuous_casting(var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp):
    h_new=[0]*var_SCZ_num # 求换热系数时的灵敏度系数
    dh=[0]*var_SCZ_num # 换热系数增量的初值
    _temp=[0]*var_SCZ_num # 计算温度初始化
    mean_temp_sensive=[0]*var_SCZ_num # 需计算温度的灵敏度系数
    errorTemp=[0]*var_SCZ_num # 表 面温度计算值与测量值之间的偏差
    initialTemp=[([var_castingTemp]*var_XNumber) for i in range(var_XNumber)]#这是一个var_XNumber*var_XNumber的列表的值为var_castingTemp
    # 整个铸坯温度场的初始温度var_XNumber=20 # 铸坯厚度方向的网格点数var_castingTemp=1530 # 铸坯的浇筑温度
    return h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp
def initialize_parameters_for_optimization(var_iter_max,var_SCZ_num):#最优化
    J=[0]*var_iter_max # 迭代过程中铸坯表面温度的测量值与计算值之间的平方差的值
    eye=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 单位矩阵的初始化
    for i in range(var_SCZ_num):
        for j in range(var_SCZ_num):
            if i==j:
                eye[i][j]=1.0
    delth=1 # 换热系数的增量
    dudh_1=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 当前换热系数得到的表面温度组成的矩阵初始化
    dudh_2=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 换热系数加入增量以后得到的表面温度组成的矩阵初始化
    dudh=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 灵敏度系数矩阵初始化
    dudh_transport=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 灵敏度系数矩阵转置初始化
    h_delt_all=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 换热系数加入增量以后得到的矩阵
    return J,eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all
def two_dimensional_mpc_solver(var_SCZ_num,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_aimTemp,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_h_initial,var_latentHeatofSolidification,var_sliceNumber,var_Z):
    h_end=[] # 要求的换热系数
    centerPoint = int(var_YNumber/2)  # 铸坯宽度截面中心点的温度
    next_temp=[([0]*var_sliceNumber) for i in range(var_sliceNumber)]# 121231拉坯方向网格划分的数量
    ########## initialize for the parameters of continuous casting start #############
    h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp=initialize_parameters_for_continuous_casting(var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp)
        ########## initialize for the parameters of continuous casting end #############
    MiddleTemp=[([var_castingTemp]*var_YNumber) for i in range(var_XNumber)]
        ########## initialize for the parameters of optimization start #############
    J,eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all=initialize_parameters_for_optimization(var_iter_max,var_SCZ_num)
        ########## initialize for the parameters of optimization end #############
        ########## equations are used to calculate the temperature of slab start #############
    tl=[0]*len(var_dis);  # 铸坯凝固时间的初值var_dis=[0.0,0.9,1.27,3.12,5.32,7.68] # 连铸二冷区各段距弯月面的距离 单位m
    for i in range(len(var_dis)):
        tl[i]=var_dis[i]/var_VcastOriginal#var_VcastOriginal=0.6/60
    time_Mold = int((tl[1]-tl[0])/var_deltTime)
    time_SCZ = int((tl[len(var_dis)-1]-tl[1])/var_deltTime)#var_deltTime=0.4 # 差分计算时间间隔
    Time_all = time_Mold+time_SCZ
    MiddleTemp_all=[([([var_castingTemp]*Time_all)]*var_YNumber) for i in range(var_XNumber)]#x,y,t三维网格初始化
    MiddleTemp_all= steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial,var_sliceNumber,var_castingTemp)
    #for stepTime in range(1,Time_all):
    #    print(MiddleTemp_all[var_XNumber-1][centerPoint][stepTime])
    slice_temp_initial =[([([var_castingTemp]*Time_all)]*var_YNumber) for i in range(var_XNumber)]
    slice_temp_all =[([([var_castingTemp]*Time_all)]*var_YNumber) for i in range(var_XNumber)]
    deltSliceNumber=0
    sliceNum=0
    delt_Slice = var_Z/var_sliceNumber    
    sliceRunningTimeAll = int((delt_Slice/var_VcastOriginal)/var_deltTime) #切片间隔内运行的时间步长的次数
    while deltSliceNumber<Time_all:
        for i in range(var_XNumber):
            for j in range(var_YNumber):
                slice_temp_initial[i][j][sliceNum]=MiddleTemp_all[i][j][deltSliceNumber]
        deltSliceNumber=deltSliceNumber+sliceRunningTimeAll
        sliceNum=sliceNum+1
        if sliceNum==var_sliceNumber:
            break
        ########## equations are used to calculate the temperature of slab end #############
   
    #补齐数据残缺
    var_runningTime=20
    var_controlTime=4
    middle_temp = [([1530]*var_XNumber) for i in range(var_XNumber)]
    var_ZNumber = 500
    #补齐完毕
    runningNumber=int(var_runningTime/var_controlTime)
    delete_create_slice_number=0
    for timeStep in range(runningNumber):
        moveSliceNumber=int(((timeStep+1)*var_controlTime)/(delt_Slice/var_VcastOriginal))
        for iter_num in range(var_iter_max):
            mean_temp,slice_temp_all=two_densonal_lifecycle_for_mpc(var_deltTime, middle_temp, var_XNumber, var_YNumber, var_X, var_Y, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_ZNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_initial,moveSliceNumber,timeStep)
            #print(mean_temp)    
                ####  model start ####
            for i in range(var_SCZ_num):
                J[iter_num]=J[iter_num]+(mean_temp[i]-var_aimTemp[i])*(mean_temp[i]-var_aimTemp[i])
                    #### model end #### where mean_temp is calculated by equation
            #print(J[iter_num])
            if J[iter_num]<var_eps:
                break
                    #### optimation start ############
            for i in range(var_SCZ_num):
                for j in range(var_SCZ_num):
                    dudh_1[i][j]=mean_temp[j]
            for i in range(var_SCZ_num):
                for j in range(var_SCZ_num):
                    if i==j:
                        h_delt_all[i][j]=var_h_initial[j]+delth
                    else:
                        h_delt_all[i][j]=var_h_initial[j]
            for i in range(4):     
                for j in range(4):
                    h_new[j]=h_delt_all[i][j]
                        ########## equations are used to calculate the temperature of slab start for sentive matrix start #############
                mean_temp_sensive,slice_temp_all=two_densonal_lifecycle_for_mpc(var_deltTime, middle_temp, var_XNumber, var_YNumber, var_X, var_Y, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_ZNumber,var_Z,var_dis,h_new,var_castingTemp,delete_create_slice_number,slice_temp_initial,moveSliceNumber,timeStep)
                #print(mean_temp_sensive)
                for k in range(var_SCZ_num):
                    dudh_2[i][k]=mean_temp_sensive[k]
            for i in range(var_SCZ_num):
                for j in range(var_SCZ_num):
                    dudh[i][j]=(dudh_1[i][j]-dudh_2[i][j])/delth
            #print(dudh)
            for i in range(var_SCZ_num):
                errorTemp[i]=(mean_temp[i]-var_aimTemp[i])
            dudh_transport=nlg.inv(eye+np.transpose(dudh)*dudh)*dudh
            for i in range(var_SCZ_num):
                dh[i]=0
                for j in range(var_SCZ_num):
                    dh[i]=dh[i]+dudh_transport[i][j]*errorTemp[j]
            for r in range(var_SCZ_num):
                var_h_initial[r]=var_h_initial[r]+dh[r]
                ########## equations are used to calculate the temperature of slab start #############
        delete_create_slice_number=moveSliceNumber
                ########## equations are used to calculate the temperature of slab end #############
    h_end = var_h_initial
    #### optimation end ############
    return h_end
