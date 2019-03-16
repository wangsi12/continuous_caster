from one_densonal_diff_lifecycle_for_mpc import one_densonal_lifecycle_for_mpc
from one_densonal_diff_lifecycle_for_mpc import steady_temp_cal
import numpy.linalg as nlg
import numpy as np
import math
def initialize_parameters_for_continuous_casting(var_ZNumber,var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp):
    h_new=[0]*var_SCZ_num # 求换热系数时的灵敏度系数
    dh=[0]*var_SCZ_num # 换热系数增量的初值
    mean_temp=[0]*var_SCZ_num # 计算温度初始化
    mean_temp_sensive=[0]*var_SCZ_num # 需计算温度的灵敏度系数
    errorTemp=[0]*var_SCZ_num # 表面温度计算值与测量值之间的偏差
    initialTemp=[([var_castingTemp]*var_ZNumber) for i in range(var_XNumber)] # 整个铸坯温度场的初始温度
    return h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp
def initialize_parameters_for_optimization(var_iter_max,var_SCZ_num):
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
    return eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all
def one_dimensional_mpc_solver(var_controlTime,var_runningTime,var_ZNumber,var_SCZ_num,var_XNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_aimTemp,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_water_initial,var_latentHeatofSolidification,var_sliceNumber,var_Z,var_VcastNew,var_alfa):
    var_h_initial=[0]*var_SCZ_num # 求换热系数
    next_temp=[([0]*var_ZNumber) for i in range(var_XNumber)]
    ########## initialize for the parameters of continuous casting start #############
    h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp=initialize_parameters_for_continuous_casting(var_ZNumber,var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp)
        ########## initialize for the parameters of continuous casting end #############
    middle_temp=[var_castingTemp]*var_ZNumber
        ########## initialize for the parameters of optimization start #############
    eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all=initialize_parameters_for_optimization(var_iter_max,var_SCZ_num)
        ########## initialize for the parameters of optimization end ############
    lengthofEachSCZ = [0]*(var_SCZ_num)
    for j in range(var_SCZ_num):
        lengthofEachSCZ[j] = var_dis[j + 2] - var_dis[j + 1]
    for k in range(var_SCZ_num):
        var_h_initial[k] = 1570 * math.pow((var_water_initial[k] / 60) / (lengthofEachSCZ[k] * var_Y * 4), 0.55) * (1 - 0.0075 * var_temperatureWater) / var_alfa[k]
        ########## equations are used to calculate the temperature of slab start #############
    tl=[0]*len(var_dis);  # 铸坯凝固时间的初值
    for i in range(len(var_dis)):
        tl[i]=var_dis[i]/var_VcastOriginal
    time_Mold = int((tl[1]-tl[0])/var_deltTime)
    time_SCZ = int((tl[len(var_dis)-1]-tl[1])/var_deltTime)
    Time_all = time_Mold+time_SCZ
    MiddleTemp_all=[[0]*(Time_all) for i in range(var_XNumber)]
    MiddleTemp_all= steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,middle_temp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial)
    slice_temp_initial =[([var_castingTemp]*var_sliceNumber) for i in range(var_XNumber)]
    slice_temp_all =[([var_castingTemp]*var_sliceNumber) for i in range(var_XNumber)]
    deltSliceNumber=0
    sliceNum=0
    delt_Slice = var_Z/var_sliceNumber    
    sliceRunningTimeAll = int((delt_Slice/var_VcastOriginal)/var_deltTime) #切片间隔内运行的时间步长的次数
    while deltSliceNumber<Time_all:
        for i in range(var_XNumber):
            slice_temp_initial[i][sliceNum]=MiddleTemp_all[i][deltSliceNumber]
        deltSliceNumber=deltSliceNumber+sliceRunningTimeAll
        sliceNum=sliceNum+1
        if sliceNum==var_sliceNumber:
            break
        ########## equations are used to calculate the temperature of slab end #############
    runningNumber=int(var_runningTime/var_controlTime)
    delete_create_slice_number=0
    water_flowerRate_all=[([0]*var_SCZ_num) for i in range(runningNumber)] # 所有的待定参数
    for timeStep in range(runningNumber):
        moveSliceNumber=int(((timeStep+1)*var_controlTime)/(delt_Slice/var_VcastNew))
        for iter_num in range(var_iter_max):
            J=[0]*var_iter_max # 迭代过程中铸坯表面温度的测量值与计算值之间的平方差的值
            mean_temp,slice_temp_all,delete_create_slice_number=one_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_X, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_ZNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_initial,var_VcastNew,moveSliceNumber,timeStep)
            #print(mean_temp)
                ####  model start ####
            for i in range(var_SCZ_num):
                J[iter_num]=J[iter_num]+(mean_temp[i]-var_aimTemp[i])*(mean_temp[i]-var_aimTemp[i])
                    #### model end #### where mean_temp is calculated by equation
            print(J[iter_num])
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
            for i in range(var_SCZ_num):     
                for j in range(var_SCZ_num):
                    h_new[j]=h_delt_all[i][j]
                        ########## equations are used to calculate the temperature of slab start for sentive matrix start #############
                mean_temp_sensive,slice_temp_all,delete_create_slice_number=one_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_X, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_ZNumber,var_Z,var_dis,h_new,var_castingTemp,delete_create_slice_number,slice_temp_initial,var_VcastNew,moveSliceNumber,timeStep)
                #print(mean_temp_sensive)
                for q in range(var_SCZ_num):
                    dudh_2[i][q]=mean_temp_sensive[q]
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
        for k in range(var_sliceNumber):
            for i in range(var_XNumber):
                slice_temp_initial[i][k]=slice_temp_all[i][k]
        delete_create_slice_number=moveSliceNumber
                ########## equations are used to calculate the temperature of slab end #############
        for k in range(no_MeasuredPoint,var_SCZ_num):
            water_flowerRate_all[timeStep, k] = math.pow((var_h_initial[k] * var_alfa[k]) / (1570 * (1 - 0.0075 * var_temperatureWater)), 1 / 0.55) * lengthofEachSCZ[k] * var_Y * 4 * 60
    #### optimation end ############
    return h_end
######### 需要输入的参数 start #################
