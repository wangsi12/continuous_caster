import math
def one_densonal_diff(h, var_deltTime, middle_temp, var_XNumber, var_X, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification):
    next_temp = [0]*var_XNumber
    deltX = var_X/var_XNumber  # 铸坯在厚度方向的空间间隔
    BOLTZMAN = 0.000000056684  # 玻尔兹曼常数
    EMISSIVITY = 0.8  # 辐射系数
    #### equation strat #####
    for i in range(var_XNumber):
        if middle_temp[i] >= var_liqTemp:
            rou = var_rouL
            specificHeat = var_specificHeatL
            Tconductivity = var_TconductivityL
        if middle_temp[i] <= var_SodTemp:
            rou = var_rouS
            specificHeat = var_specificHeatS
            Tconductivity = var_TconductivityS
        if (var_SodTemp < middle_temp[i]) & (var_liqTemp > middle_temp[i]):
            rou = (var_rouS-var_rouL)*(var_liqTemp -
                                       middle_temp[i])/(var_liqTemp-var_SodTemp)+var_rouL
            Tconductivity = (var_TconductivityS)*(var_liqTemp-middle_temp[i])/(var_liqTemp-var_SodTemp)+var_m*(
                1-(var_liqTemp-middle_temp[i])/(var_liqTemp-var_SodTemp))*var_TconductivityL
            specificHeat = (var_specificHeatS-var_specificHeatL)*(var_liqTemp-middle_temp[i])/(
                var_liqTemp-var_SodTemp)+var_specificHeatL+var_latentHeatofSolidification/(var_liqTemp-var_SodTemp)
        a1 = Tconductivity/(rou*specificHeat)
        a = (Tconductivity*var_deltTime)/(rou*specificHeat*deltX*deltX)
        if i == 0:
            next_temp[i] = (1-2*var_deltTime*h/(rou*specificHeat*deltX)-2*var_deltTime*a1/(deltX*deltX))*middle_temp[i] + \
                2*var_deltTime*a1 * \
                middle_temp[i+1]/(deltX*deltX)+2*var_deltTime * \
                h*var_temperatureWater/(rou*specificHeat*deltX)
        if i == var_XNumber-1:
            next_temp[i] = (1-2*var_deltTime*h/(rou*specificHeat*deltX)-2*var_deltTime*a1/(deltX*deltX))*middle_temp[i] + \
                2*var_deltTime*a1 * \
                middle_temp[i-1]/(deltX*deltX)+2*var_deltTime * \
                h*var_temperatureWater/(rou*specificHeat*deltX)
        if 0 < i < var_XNumber-1:
            next_temp[i] = a*middle_temp[i+1]+a * \
                middle_temp[i-1]+(1-2*a)*middle_temp[i]
    #### equation end #####
    return next_temp
def select_h_in_different_section(k,var_ZNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastOriginal):
    delt_z=var_Z/var_ZNumber
    section_location_now=delt_z*k
    h=1.0
    if var_dis[0]<=section_location_now<=var_dis[1]:
        tTime = section_location_now/var_VcastOriginal;
        h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
        sectionNumber[0]=sectionNumber[0]+1
    if var_dis[1]<section_location_now<=var_dis[2]:
        h=var_h_initial[0]
        sectionNumber[1]=sectionNumber[1]+1          
    if var_dis[2]<section_location_now<=var_dis[3]:
        h=var_h_initial[1]
        sectionNumber[2]=sectionNumber[2]+1
    if var_dis[3]<section_location_now<=var_dis[4]:
        h=var_h_initial[2]
        sectionNumber[3]=sectionNumber[3]+1
    if var_dis[4]<section_location_now<=var_dis[5]:
        h=var_h_initial[3]
        sectionNumber[4]=sectionNumber[4]+1
    if var_dis[5]<section_location_now<=var_dis[6]:
        h=var_h_initial[4]
        sectionNumber[5]=sectionNumber[5]+1
    if var_dis[6]<section_location_now<=var_dis[7]:
        h=var_h_initial[5]
        sectionNumber[6]=sectionNumber[6]+1
    if var_dis[7]<section_location_now<=var_dis[8]:
        h=var_h_initial[6]
        sectionNumber[7]=sectionNumber[7]+1
    if var_dis[8]<section_location_now:
        h=var_h_initial[7]
        sectionNumber[8]=sectionNumber[8]+1
    return h,sectionNumber
def steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial):
    NextTemp=[0]*(var_XNumber)
    MiddleTemp_all=[[0]*(Time_all) for i in range(var_XNumber)]
    for stepTime in range(1,Time_all):
        if stepTime <= time_Mold:
            tTime = var_deltTime*stepTime;
            h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
            NextTemp=one_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification)
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
            NextTemp=one_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification)
        for i in range (var_XNumber):
            MiddleTemp[i]=NextTemp[i]
            MiddleTemp_all[i][stepTime]=NextTemp[i]
    return MiddleTemp_all
def one_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_X, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_sliceNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_all,var_VcastNew,moveSliceNumber,stepTime_for_running):
    delt_Slice = var_Z/var_sliceNumber
    slice_temp_all_middle =[([var_castingTemp]*var_sliceNumber) for i in range(var_XNumber)]
    slice_temp_new =[([var_castingTemp]*var_sliceNumber) for i in range(var_XNumber)]
    for k in range(var_sliceNumber):
        for j in range(var_XNumber):
            slice_temp_new[j][k]=slice_temp_all[j][k]
    next_temp = [0]*var_XNumber
    middle_temp = [0]*var_XNumber
    mean_temp=[0]*(len(var_dis)-2) # 平均温度的初始化
    runningTimeinControltime = int(var_controlTime/var_deltTime)  #控制周期内的迭代次数
    sliceRunningTimeAll_New = int((delt_Slice/var_VcastNew)/var_deltTime) #切片间隔内运行的时间步长的次数
    delt_delete_create_slice_number=moveSliceNumber-delete_create_slice_number
    print(var_h_initial)
    for deltSliceNum in range(delt_delete_create_slice_number):
        sectionNumber=[0]*(len(var_dis)-1) # 连铸二冷区各段切片数量
        for k in range(var_sliceNumber):
            h,sectionNumber=select_h_in_different_section(k,var_sliceNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastNew)
            for i in range(var_XNumber):
                middle_temp[i]=slice_temp_new[i][k]
            for deltTime_Number in range(sliceRunningTimeAll_New):
                next_temp=one_densonal_diff(h,var_deltTime,middle_temp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification)
                for q in range(var_XNumber):
                    middle_temp[q]=next_temp[q]
            for i in range(var_XNumber):
                slice_temp_new[i][k]=middle_temp[i]
                slice_temp_all_middle[i][k]=middle_temp[i]
        for k in range(var_sliceNumber):
            if k==0:
                for i in range(var_XNumber):
                    slice_temp_new[i][0]=var_castingTemp
            else:
                for j in range(var_XNumber):
                    slice_temp_new[j][k]=slice_temp_all_middle[j][k-1]
        for k in range(var_sliceNumber):
            for i in range(var_XNumber):
                slice_temp_all_middle[i][k]=slice_temp_new[i][k]
    MoldNumber = sectionNumber[0]
    sectionNumber_one = sectionNumber[0]+sectionNumber[1]
    sectionNumber_two = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]
    sectionNumber_three = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]
    sectionNumber_four = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]
    sectionNumber_five = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]+sectionNumber[5]
    sectionNumber_six = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]+sectionNumber[5]+sectionNumber[6]
    sectionNumber_seven = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]+sectionNumber[5]+sectionNumber[6]+sectionNumber[7]
    sectionNumber_eight = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]+sectionNumber[5]+sectionNumber[6]+sectionNumber[7]+sectionNumber[8]
    for i in range(MoldNumber,sectionNumber_one-1):
        mean_temp[0] = mean_temp[0] + slice_temp_all_middle[0][i]
    mean_temp[0] = mean_temp[0]/sectionNumber[1]
    for i in range(sectionNumber_one,sectionNumber_two-1):
        mean_temp[1] = mean_temp[1] + slice_temp_all_middle[0][i]
    mean_temp[1] = mean_temp[1]/sectionNumber[2]
    for i in range(sectionNumber_two,sectionNumber_three-1):
        mean_temp[2] = mean_temp[2] + slice_temp_all_middle[0][i]
    mean_temp[2] = mean_temp[2]/sectionNumber[3]
    for i in range(sectionNumber_three,sectionNumber_four-1):
        mean_temp[3] = mean_temp[3] + slice_temp_all_middle[0][i]
    mean_temp[3] = mean_temp[3]/sectionNumber[4]
    for i in range(sectionNumber_four,sectionNumber_five-1):
        mean_temp[4] = mean_temp[4] + slice_temp_all_middle[0][i]
    mean_temp[4] = mean_temp[4]/sectionNumber[5]  
    for i in range(sectionNumber_five,sectionNumber_six-1):
        mean_temp[5] = mean_temp[5] + slice_temp_all_middle[0][i]
    mean_temp[5] = mean_temp[5]/sectionNumber[6]   
    for i in range(sectionNumber_six,sectionNumber_seven-1):
        mean_temp[6] = mean_temp[6] + slice_temp_all_middle[0][i]
    mean_temp[6] = mean_temp[6]/sectionNumber[7] 
    for i in range(sectionNumber_seven,sectionNumber_eight-1):
        mean_temp[7] = mean_temp[7] + slice_temp_all_middle[0][i]
    mean_temp[7] = mean_temp[7]/sectionNumber[8]      
    return mean_temp,slice_temp_all_middle,delete_create_slice_number

    

