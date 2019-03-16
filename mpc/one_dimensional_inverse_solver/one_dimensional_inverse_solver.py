from one_dinensional_diff_solver import one_densonal_diff
import numpy.linalg as nlg
import numpy as np
import math
 
def initialize_parameters_for_continuous_casting(var_SCZ_num,var_numberofMeasuredPoint,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp,var_disMeasure):
    no_MeasuredPoint=var_SCZ_num-var_numberofMeasuredPoint # 为测量点的个数
    h_new=[0]*var_SCZ_num # 求换热系数时的灵敏度系数
    dh=[0]*var_numberofMeasuredPoint # 换热系数增量的初值
    tl=[0]*len(var_dis);  # 铸坯凝固时间的初值
    for i in range(len(var_dis)):
        tl[i]=var_dis[i]/var_VcastOriginal
    time_Mold = int((tl[1]-tl[0])/var_deltTime)
    time_SCZ = int((tl[len(var_dis)-1]-tl[1])/var_deltTime)
    Time_all = time_Mold+time_SCZ
    measuredStepPoint=[0]*var_numberofMeasuredPoint # 网格划分以后测量温度在网格中的位置
    for i in range(var_numberofMeasuredPoint):
        measuredStepPoint[i]=int((var_disMeasure[i]/var_VcastOriginal)/var_deltTime)
    calTemp=[0]*var_numberofMeasuredPoint # 计算温度初始化
    calTemp_sensive=[0]*var_numberofMeasuredPoint # 需计算温度的灵敏度系数
    errorTemp=[0]*var_numberofMeasuredPoint # 表面温度计算值与测量值之间的偏差
    initialTemp=[1]*var_XNumber # 整个铸坯温度场的初始温度
    for i in range(var_XNumber):
        initialTemp[i]=var_castingTemp
    return no_MeasuredPoint,h_new,dh,time_Mold,time_SCZ,Time_all,measuredStepPoint,calTemp,calTemp_sensive,errorTemp,initialTemp
def initialize_parameters_for_optimization(var_iter_max,var_numberofMeasuredPoint,var_SCZ_num):
    J=[0]*var_iter_max # 迭代过程中铸坯表面温度的测量值与计算值之间的平方差的值
    eye=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 单位矩阵的初始化
    for i in range(var_numberofMeasuredPoint):
        for j in range(var_numberofMeasuredPoint):
            if i==j:
                eye[i][j]=1.0
    delth=1 # 换热系数的增量
    dudh_1=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 当前换热系数得到的表面温度组成的矩阵初始化
    dudh_2=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 换热系数加入增量以后得到的表面温度组成的矩阵初始化
    dudh=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 灵敏度系数矩阵初始化
    dudh_transport=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 灵敏度系数矩阵转置初始化
    h_delt_all=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 换热系数加入增量以后得到的矩阵
    return J,eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all
def temp_cal_in_measured_point(var_numberofMeasuredPoint,var_h_initial,var_VcastOriginal,var_dis,initialTemp,Time_all,time_Mold,var_deltTime,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification,measuredStepPoint):
    calTemp=[0]*var_numberofMeasuredPoint # 计算温度初始化
    TemSur = [0]*(Time_all) # 表面温度的初值
    MiddleTemp=[0]*var_XNumber
    NextTemp=[0]*var_XNumber
    for i in range(var_XNumber):
        MiddleTemp[i]=initialTemp[i] # 当前时刻铸坯界面的温度分布
        NextTemp[i]=initialTemp[i] # 下一时刻铸坯界面的温度分布
    for stepTime in range(1,Time_all):
        if stepTime <= time_Mold:
            tTime = var_deltTime*stepTime;
            h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
            NextTemp=one_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification)
        else:
            disNow = var_dis[1]+stepTime*var_VcastOriginal*var_deltTime
            if var_dis[1]<=disNow<=var_dis[2]:
                h=var_h_initial[0]          
            if var_dis[2]<disNow<=var_dis[3]:
                h=var_h_initial[1]
            if var_dis[3]<disNow<=var_dis[4]:
                h=var_h_initial[2]
            if var_dis[4]<disNow<=var_dis[5]:
                h=var_h_initial[3]
            if var_dis[5]<disNow<=var_dis[6]:
                h=var_h_initial[4]
            if var_dis[6]<disNow<=var_dis[7]:
                h=var_h_initial[5]
            if var_dis[7]<disNow<=var_dis[8]:
                h=var_h_initial[6]
            if var_dis[8]<disNow:
                h=var_h_initial[7]
            NextTemp=one_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification)
        for i in range (var_XNumber):
            MiddleTemp[i]=NextTemp[i]
        TemSur[stepTime-1]=MiddleTemp[0]
        for j in range(var_numberofMeasuredPoint):
            calTemp[j]=TemSur[measuredStepPoint[j]]
    return calTemp

def one_dimensional_inverse_solver(var_XNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_measuredTemp,var_disMeasure,var_numberofMeasuredPoint,var_SCZ_num,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_alfa_initial,latentHeatofSolidification,var_WaterFlowRate,var_NumberofSCZ):
    alfa_get=[0]*var_SCZ_num # 待定参数的初值
    var_h_initial=[0]*var_SCZ_num # 求换热系数
    ########## initialize for the parameters of continuous casting start #############
    no_MeasuredPoint,h_new,dh,time_Mold,time_SCZ,Time_all,measuredStepPoint,calTemp,calTemp_sensive,errorTemp,initialTemp=initialize_parameters_for_continuous_casting(var_SCZ_num,var_numberofMeasuredPoint,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp,var_disMeasure)
        ########## initialize for the parameters of continuous casting end #############

        ########## initialize for the parameters of optimization start #############
    J,eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all=initialize_parameters_for_optimization(var_iter_max,var_numberofMeasuredPoint,var_SCZ_num)
        ########## initialize for the parameters of optimization end #############
    lengthofEachSCZ = [0]*(var_NumberofSCZ)
    for j in range(var_NumberofSCZ):
        lengthofEachSCZ[j] = var_dis[j + 2] - var_dis[j + 1]
    for k in range(no_MeasuredPoint,var_SCZ_num):
        var_h_initial[k] = 1570 * math.pow((var_WaterFlowRate[k] / 60) / (lengthofEachSCZ[k] * var_Y * 4), 0.55) * (1 - 0.0075 * var_temperatureWater) / var_alfa_initial[k]
        ########## equations are used to calculate the temperature of slab start #############
    calTemp=temp_cal_in_measured_point(var_numberofMeasuredPoint,var_h_initial,var_VcastOriginal,var_dis,initialTemp,Time_all,time_Mold,var_deltTime,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification,measuredStepPoint)
        ########## equations are used to calculate the temperature of slab end #############
    for iter_num in range(var_iter_max):
                ####  model start ####
        for i in range(var_numberofMeasuredPoint):
            J[iter_num]=J[iter_num]+(calTemp[i]-var_measuredTemp[i])*(calTemp[i]-var_measuredTemp[i])
                #### model end #### where calTemp is calculated by equation
        print(J[iter_num])
        if J[iter_num]<var_eps:
            break
            #J[iter_num]=0
                #### optimation start ############
        for i in range(var_numberofMeasuredPoint):
            for j in range(var_numberofMeasuredPoint):
                dudh_1[i][j]=calTemp[j]
        for i in range(var_SCZ_num):
            for j in range(var_SCZ_num):
                if i==j:
                    h_delt_all[i][j]=var_h_initial[j]+delth
                else:
                    h_delt_all[i][j]=var_h_initial[j]
        for i in range(no_MeasuredPoint):
            h_delt_all[i][i]=h_delt_all[i][i]-delth
        for i in range(no_MeasuredPoint,var_SCZ_num):     
            for j in range(var_SCZ_num):
                h_new[j]=h_delt_all[i][j]
                    ########## equations are used to calculate the temperature of slab start for sentive matrix start #############
                calTemp_sensive=temp_cal_in_measured_point(var_numberofMeasuredPoint,h_new,var_VcastOriginal,var_dis,initialTemp,Time_all,time_Mold,var_deltTime,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification,measuredStepPoint)
                    ########## equations are used to calculate the temperature of slab start for sentive matrix end #############
            for k in range(var_numberofMeasuredPoint):
                dudh_2[i-no_MeasuredPoint][k]=calTemp_sensive[k]
        for i in range(var_numberofMeasuredPoint):
            for j in range(var_numberofMeasuredPoint):
                dudh[i][j]=(dudh_1[i][j]-dudh_2[i][j])/delth
        for i in range(var_numberofMeasuredPoint):
            errorTemp[i]=(calTemp[i]-var_measuredTemp[i])
        dudh_transport=nlg.inv(eye+np.transpose(dudh)*dudh)*dudh
        for i in range(var_numberofMeasuredPoint):
            dh[i]=0
            for j in range(var_numberofMeasuredPoint):
                dh[i]=dh[i]+dudh_transport[i][j]*errorTemp[j]
        for r in range(no_MeasuredPoint,var_SCZ_num):
            var_h_initial[r]=var_h_initial[r]+dh[r-no_MeasuredPoint]
            ########## equations are used to calculate the temperature of slab start #############
        calTemp=temp_cal_in_measured_point(var_numberofMeasuredPoint,var_h_initial,var_VcastOriginal,var_dis,initialTemp,Time_all,time_Mold,var_deltTime,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,latentHeatofSolidification,measuredStepPoint)
            ########## equations are used to calculate the temperature of slab end #############
    for k in range(no_MeasuredPoint,var_SCZ_num):
        alfa_get[k] = 1570 * math.pow((var_WaterFlowRate[k] / 60) / (lengthofEachSCZ[k] * Y * 4), 0.55) * (1 - 0.0075 * var_temperatureWater) / var_h_initial[k]
    #### optimation end ############
    return alfa_get
