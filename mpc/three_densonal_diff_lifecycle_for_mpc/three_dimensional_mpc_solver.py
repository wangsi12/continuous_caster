from three_densonal_diff_lifecycle_for_mpc import three_densonal_lifecycle_for_mpc
import numpy.linalg as nlg
import numpy as np
def initialize_parameters_for_continuous_casting(var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp,var_sliceNumber):
    h_new=[0]*var_SCZ_num # 求换热系数时的灵敏度系数
    dh=[0]*var_SCZ_num # 换热系数增量的初值
    mean_temp=[0]*var_SCZ_num # 计算温度初始化
    mean_temp_sensive=[0]*var_SCZ_num # 需计算温度的灵敏度系数
    errorTemp=[0]*var_SCZ_num # 表 面温度计算值与测量值之间的偏差
    initialTemp=[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]# 整个铸坯温度场的初始温度
    return h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp
def initialize_parameters_for_optimization(var_iter_max,var_SCZ_num):
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
def three_dimensional_mpc_solver(var_SCZ_num,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_aimTemp,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_h_initial,var_latentHeatofSolidification,var_sliceNumber,var_Z,Var_runningSteadyTime,var_VcastNew):
    h_end=[] # 要求的换热系数
    next_temp=[([0]*var_ZNumber) for i in range(var_XNumber)]
    ########## initialize for the parameters of continuous casting start #############
    h_new,dh,mean_temp,mean_temp_sensive,errorTemp,initialTemp=initialize_parameters_for_continuous_casting(var_SCZ_num,var_dis,var_VcastOriginal,var_deltTime,var_XNumber,var_castingTemp,var_sliceNumber)
        ########## initialize for the parameters of continuous casting end #############
    middle_temp=initialTemp
        ########## initialize for the parameters of optimization start #############
    J,eye,delth,dudh_1,dudh_2,dudh,dudh_transport,h_delt_all=initialize_parameters_for_optimization(var_iter_max,var_SCZ_num)
        ########## initialize for the parameters of optimization end #############
        ########## equations are used to calculate the temperature of slab start #############
    steadyRunningNumber=int(Var_runningSteadyTime/var_controlTime)
    slice_temp_all_middle =[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    slice_temp_initial =[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    delete_create_slice_number=0
    delt_Slice = var_Z/var_sliceNumber
    for steadyStepTime in range(steadyRunningNumber):
        moveSliceNumber=int(((steadyStepTime+1)*var_controlTime)/(delt_Slice/var_VcastOriginal))
        mean_temp,slice_temp_all_middle=three_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_ZNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_initial,moveSliceNumber,steadyStepTime)
        delete_create_slice_number=moveSliceNumber
        print(mean_temp)
        for k in range(var_sliceNumber):
            for i in range(var_YNumber):
                for j in range(var_XNumber):
                    slice_temp_initial[i][j][k]=slice_temp_all_middle[i][j][k]
        ########## equations are used to calculate the temperature of slab end #############
    runningNumber=int(var_runningTime/var_controlTime)
    delete_create_slice_number=0
    for timeStep in range(runningNumber):
        moveSliceNumber=int(((timeStep+1)*var_controlTime)/(delt_Slice/var_VcastNew))
        for iter_num in range(var_iter_max):
            mean_temp,slice_temp_all_middle=three_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_YNumber, var_X, var_Y, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastNew,var_controlTime,var_ZNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_initial,moveSliceNumber,timeStep)
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
                mean_temp_sensive,slice_temp_all_middle=three_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_YNumber, var_X, var_Y, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastNew,var_controlTime,var_ZNumber,var_Z,var_dis,h_new,var_castingTemp,delete_create_slice_number,slice_temp_initial,moveSliceNumber,timeStep)
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
######### 需要输入的参数 start #################
