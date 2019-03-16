import math
def two_densonal_diff(h, var_deltTime, middle_temp, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification):
    deltX = var_X/var_XNumber  # 铸坯在厚度方向的空间间隔
    deltY = var_Y/var_YNumber  # 铸坯在宽度方向的空间间隔
    BOLTZMAN = 0.000000056684  # 玻尔兹曼常数
    EMISSIVITY = 0.8  # 辐射系数
    next_temp = [[0]*var_YNumber for i in range(var_XNumber)]
    #### equation strat #####
    for i in range(var_XNumber):
        for j in range(var_YNumber):
            if middle_temp[i][j] >= var_liqTemp:
                rou = var_rouL
                specificHeat = var_specificHeatL
                Tconductivity = var_TconductivityL
            if middle_temp[i][j] <= var_SodTemp:
                rou = var_rouS
                specificHeat = var_specificHeatS
                Tconductivity = var_TconductivityS
            if (var_SodTemp < middle_temp[i][j]) & (var_liqTemp > middle_temp[i][j]):
                rou = (var_rouS-var_rouL)*(var_liqTemp -
                                        middle_temp[i][j])/(var_liqTemp-var_SodTemp)+var_rouL
                Tconductivity = (var_TconductivityS)*(var_liqTemp-middle_temp[i][j])/(var_liqTemp-var_SodTemp)+var_m*(
                    1-(var_liqTemp-middle_temp[i][j])/(var_liqTemp-var_SodTemp))*var_TconductivityL
                specificHeat = (var_specificHeatS-var_specificHeatL)*(var_liqTemp-middle_temp[i][j])/(
                    var_liqTemp-var_SodTemp)+var_specificHeatL+var_latentHeatofSolidification/(var_liqTemp-var_SodTemp)
            # a1 = Tconductivity/(rou*specificHeat)
            a = (Tconductivity*var_deltTime)/(rou*specificHeat*deltX*deltX)
            b = (Tconductivity*var_deltTime)/(rou*specificHeat*deltY*deltY)
            ##### 四个边的温度 start ######
            if i==0 and (j!=0 and j!=var_YNumber-1):    # 情况1
                next_temp[i][j] = middle_temp[i][j]+2*a*(middle_temp[i+1][j]-middle_temp[i][j]) + b*(middle_temp[i][j+1]-2*middle_temp[i][j]+middle_temp[i][j-1])
            if j==0 and (i!=var_XNumber-1 and i!=0):  # 情况2
                next_temp[i][j] = middle_temp[i][j]+2*a*(middle_temp[i][j+1]-middle_temp[i][j]) + b*(middle_temp[i+1][j]-2*middle_temp[i][j]+middle_temp[i-1][j])
            if i==var_XNumber-1 and (j!=var_YNumber-1 and j!=0):   # 情况3
                h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp[i][j]*middle_temp[i][j]+var_temperatureWater*var_temperatureWater)*(middle_temp[i][j]+var_temperatureWater)
                next_temp[i][j]=middle_temp[i][j]+2*a*(middle_temp[i-1][j]-middle_temp[i][j])+b*(middle_temp[i][j+1]-2*middle_temp[i][j]+middle_temp[i][j-1])-(2*h_int*deltX*a*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)
            if j==var_YNumber-1 and (i!=var_XNumber-1 and i!=0):     # 情况4
                h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp[i][j]*middle_temp[i][j]+var_temperatureWater*var_temperatureWater)*(middle_temp[i][j]+var_temperatureWater)
                next_temp[i][j]=middle_temp[i][j]+2*b*(middle_temp[i][j-1]-middle_temp[i][j])+a*(middle_temp[i+1][j]-2*middle_temp[i][j]+middle_temp[i-1][j])-(2*h_int*deltY*b*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)
            ##### 四个边的温度 end ######
            ##### 四个角的温度 start ####
            if i==0 and j==0:           # 情况5
                next_temp[i][j]=middle_temp[i][j]+2*a*(middle_temp[i+1][j]-middle_temp[i][j])+2*b*(middle_temp[i][j+1]-middle_temp[i][j])
            if i==0 and j==var_YNumber-1:      # 情况6
                h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp[i][j]*middle_temp[i][j]+var_temperatureWater*var_temperatureWater)*(middle_temp[i][j]+var_temperatureWater)
                next_temp[i][j]=middle_temp[i][j]+2*a*(middle_temp[i+1][j]-middle_temp[i][j])+2*b*(middle_temp[i][j-1]-middle_temp[i][j])-(2*h_int*deltY*b*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)
            if i==var_XNumber-1 and j==0:          # 情况7
                h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp[i][j]*middle_temp[i][j]+var_temperatureWater*var_temperatureWater)*(middle_temp[i][j]+var_temperatureWater)
                next_temp[i][j]=middle_temp[i][j]+2*a*(middle_temp[i-1][j]-middle_temp[i][j])+2*b*(middle_temp[i][j+1]-middle_temp[i][j])-(2*h_int*deltX*a*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)
            if i==var_XNumber-1 and j==var_YNumber-1:         # 情况8
                h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp[i][j]*middle_temp[i][j]+var_temperatureWater*var_temperatureWater)*(middle_temp[i][j]+var_temperatureWater)
                next_temp[i][j]=middle_temp[i][j]+2*a*(middle_temp[i-1][j]-middle_temp[i][j])+2*b*(middle_temp[i][j-1]-middle_temp[i][j])-(2*h_int*deltY*b*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)-(2*h_int*deltX*a*(middle_temp[i][j]-var_temperatureWater)/Tconductivity)
            ##### 四个角的温度 end ####
            ##### 内部温度 start ####
            if  (i!=0 and i!=var_XNumber-1) and (j!=0 and j!=var_YNumber-1):          # 情况9
                next_temp[i][j]=middle_temp[i][j]+a*(middle_temp[i+1][j]-2*middle_temp[i][j]+middle_temp[i-1][j])+b*(middle_temp[i][j+1]-2*middle_temp[i][j]+middle_temp[i][j-1])
            ##### 内部温度 end ####
    #### equation end #####
    print("一次开始")
    print( next_temp)
    print("一次结束")
    return next_temp
def select_h_in_different_section(k,var_ZNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastOriginal,deltTime_Number,var_deltTime):
    delt_z=var_Z/var_ZNumber
    section_location_now=delt_z*k
    h=1.0
    if var_dis[0]<=section_location_now<=var_dis[1]:
        tTime = section_location_now/var_VcastOriginal
        h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
        if deltTime_Number==0:
            sectionNumber[0]=sectionNumber[0]+1
    if var_dis[1]<section_location_now<=var_dis[2]:
        h=var_h_initial[0]
        if deltTime_Number==0:
            sectionNumber[1]=sectionNumber[1]+1          
    if var_dis[2]<section_location_now<=var_dis[3]:
        h=var_h_initial[1]
        if deltTime_Number==0:
            sectionNumber[2]=sectionNumber[2]+1
    if var_dis[3]<section_location_now<=var_dis[4]:
        h=var_h_initial[2]
        if deltTime_Number==0:
            sectionNumber[3]=sectionNumber[3]+1
    if var_dis[4]<section_location_now<=var_dis[5]:
        h=var_h_initial[3]
        if deltTime_Number==0:
            sectionNumber[4]=sectionNumber[4]+1
    return h,sectionNumber
def steady_temp_cal(var_dis,var_VcastOriginal,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,Time_all,time_Mold,var_h_initial,var_sliceNumber,var_castingTemp):
    NextTemp=[([0]*var_YNumber) for i in range(var_XNumber)]
    #MiddleTemp_all =[([([var_castingTemp]*Time_all)]*var_YNumber) for i in range(var_XNumber)]
    
    MiddleTemp_all = [0] * var_XNumber
    for i in range(var_XNumber):
        MiddleTemp_all[i] = [0] * var_YNumber
    for i in range(var_XNumber):
        for j in range(var_YNumber):
            MiddleTemp_all[i][j] = [0] * Time_all

    for stepTime in range(Time_all):
        if stepTime <= time_Mold:
            tTime = var_deltTime*(stepTime+1);
            h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
            NextTemp=two_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification)
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
            NextTemp=two_densonal_diff(h,var_deltTime,MiddleTemp,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification)
        for i in range (var_XNumber):
            for j in range(var_YNumber):
                MiddleTemp[i][j]=NextTemp[i][j]
                MiddleTemp_all[i][j][stepTime]=NextTemp[i][j]
    return MiddleTemp_all
def two_densonal_lifecycle_for_mpc(var_deltTime, middle_temp, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_sliceNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_all,moveSliceNumber,stepTime_for_running):
    delt_Slice = var_Z/var_sliceNumber
    slice_temp_new =[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    slice_temp_all_middle=[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    for k in range(var_sliceNumber):
        for i in range(var_XNumber):
            for j in range(var_YNumber):
                slice_temp_new[i][j][k]=slice_temp_all[i][j][k]
    next_temp = [([0]*var_YNumber) for i in range(var_XNumber)]
    mean_temp=[0]*(len(var_dis)-2) # 平均温度的初始化
    runningTimeinControltime = int(var_controlTime/var_deltTime)  #控制周期内的迭代次数
    sliceRunningTimeAll = int((delt_Slice/var_VcastOriginal)/var_deltTime) #切片间隔内运行的时间步长的次数
    centerPoint = int(var_YNumber/2)  # 铸坯宽度截面中心点的温度
    delt_delete_create_slice_number=moveSliceNumber-delete_create_slice_number
    sectionNumber=[0]*(len(var_dis)-1) # 连铸二冷区各段切片数量
    #print(moveSliceNumber,delete_create_slice_number)
    for deltSliceNum in range(delt_delete_create_slice_number):
        for i in range(len(var_dis)-1):
            sectionNumber[i]=0 # 连铸二冷区各段切片数量
        for k in range(var_sliceNumber):
            for i in range(var_XNumber):
                for j in range(var_YNumber):
                    middle_temp[i][j]=slice_temp_new[i][j][k]
            for deltTime_Number in range(sliceRunningTimeAll):
                h,sectionNumber=select_h_in_different_section(k,var_sliceNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastOriginal,deltTime_Number,var_deltTime)
                next_temp=two_densonal_diff(h, var_deltTime, middle_temp, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m, var_latentHeatofSolidification)
            for i in range(var_XNumber):
                for j in range(var_YNumber):
                    middle_temp[i][j]=next_temp[i][j]
            for i in range(var_XNumber):
                for j in range(var_YNumber):
                    slice_temp_all_middle[i][j][k]=middle_temp[i][j]
                    slice_temp_new[i][j][k]=middle_temp[i][j]
        for k in range(var_sliceNumber):
            if k==0:
                for i in range(var_XNumber):
                    for j in range(var_YNumber):
                        slice_temp_new[i][j][0]=var_castingTemp
            else:
                for i in range(var_XNumber):
                    for j in range(var_YNumber):
                        slice_temp_new[i][j][k]=slice_temp_all_middle[i][j][k-1]
    MoldNumber = sectionNumber[0]
    sectionNumber_one = sectionNumber[0]+sectionNumber[1]
    sectionNumber_two = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]
    sectionNumber_three = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]
    sectionNumber_four = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]
    #print(sectionNumber)
    for i in range(MoldNumber,sectionNumber_one-1):
        mean_temp[0] = mean_temp[0] + slice_temp_new[var_XNumber-1][centerPoint][i]
    mean_temp[0] = mean_temp[0]/sectionNumber[1]
    for i in range(sectionNumber_one,sectionNumber_two-1):
        mean_temp[1] = mean_temp[1] + slice_temp_new[var_XNumber-1][centerPoint][i]
    mean_temp[1] = mean_temp[1]/sectionNumber[2]
    for i in range(sectionNumber_two,sectionNumber_three-1):
        mean_temp[2] = mean_temp[2] + slice_temp_new[var_XNumber-1][centerPoint][i]
    mean_temp[2] = mean_temp[2]/sectionNumber[3]
    for i in range(sectionNumber_three,sectionNumber_four-1):
        mean_temp[3] = mean_temp[3] + slice_temp_new[var_XNumber-1][centerPoint][i]
    mean_temp[3] = mean_temp[3]/sectionNumber[4]   
    return mean_temp,slice_temp_new

    

