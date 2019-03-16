import math
def three_densonal_diff(var_deltTime, middle_temp_all, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_sliceNumber,var_Z,var_dis,var_h_initial,var_VcastOriginal):
    deltX = var_X/var_XNumber  # 铸坯在厚度方向的空间间隔
    deltY = var_Y/var_YNumber  # 铸坯在宽度方向的空间间隔
    deltZ = var_Z/var_sliceNumber # 铸坯在拉坯方向的空间间隔
    BOLTZMAN = 0.000000056684  # 玻尔兹曼常数
    EMISSIVITY = 0.8  # 辐射系数
    next_temp_all = [([([1]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    sectionNumber=[0]*(len(var_dis)-1) # 连铸二冷区各段切片数量
    #### equation strat #####
    for k in range(var_sliceNumber):
        h,sectionNumber=select_h_in_different_section(k,var_sliceNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastOriginal)
        for i in range(var_XNumber):
            for j in range(var_YNumber):
                if middle_temp_all[i][j][k] >= var_liqTemp:
                    rou = var_rouL
                    specificHeat = var_specificHeatL
                    Tconductivity = var_TconductivityL
                if middle_temp_all[i][j][k] <= var_SodTemp:
                    rou = var_rouS
                    specificHeat = var_specificHeatS
                    Tconductivity = var_TconductivityS
                if (var_SodTemp < middle_temp_all[i][j][k]) & (var_liqTemp > middle_temp_all[i][j][k]):
                    rou = (var_rouS-var_rouL)*(var_liqTemp -
                                            middle_temp_all[i][j][k])/(var_liqTemp-var_SodTemp)+var_rouL
                    Tconductivity = (var_TconductivityS)*(var_liqTemp-middle_temp_all[i][j][k])/(var_liqTemp-var_SodTemp)+var_m*(
                        1-(var_liqTemp-middle_temp_all[i][j][k])/(var_liqTemp-var_SodTemp))*var_TconductivityL
                    specificHeat = (var_specificHeatS-var_specificHeatL)*(var_liqTemp-middle_temp_all[i][j][k])/(
                        var_liqTemp-var_SodTemp)+var_specificHeatL+var_latentHeatofSolidification/(var_liqTemp-var_SodTemp)
                # a1 = Tconductivity/(rou*specificHeat)
                a = (Tconductivity*var_deltTime)/(rou*specificHeat*deltX*deltX)
                b = (Tconductivity*var_deltTime)/(rou*specificHeat*deltY*deltY)
                d = (Tconductivity*var_deltTime)/(rou*specificHeat*deltZ*deltZ)
                ####### 弯月面处（z=0）网格节点的温度 strat ################
                ##########四个边 start  ###########
                if i==0 and (j!=0 and j!=var_YNumber-1) and k==0: # 情况1 左边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and (j!=0 and j!=var_YNumber-1) and k==0: # 情况2 右边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if j==0 and (i!=0 and i!=var_XNumber-1) and k==0: # 情况3 上边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if j==var_YNumber-1 and (i!=0 and i!=var_XNumber-1) and k==0: # 情况4 下边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ##########四个边 end    ###########
                ############ 四个角 start ###############
                if i==0 and j==0 and k==0: #########  情况5 原点O1点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==0 and k==0: #########  情况6 A1点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==0 and j==var_YNumber-1 and k==0: #########  情况7 B1点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==var_YNumber-1 and k==0: #########  情况8 C1点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ############ 四个角 end ###############
                if (i!=0 and i!=var_XNumber-1) and (j!=0 and j!=var_YNumber-1) and k==0: #########  情况9 表面温度
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k+1]-middle_temp_all[i][j][k])
                ####### 弯月面处（z=0）网格节点的温度 end ################
                ####### 二冷区出口处（z=L）网格节点的温度 start ################
                if i==0 and (j!=0 and j!=var_YNumber-1) and k==var_sliceNumber-1:# 情况10 左边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and (j!=0 and j!=var_YNumber-1) and k==var_sliceNumber-1: # 情况11 右边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if j==0 and (i!=0 and i!=var_XNumber-1) and k==var_sliceNumber-1: # 情况12 上边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if j==var_YNumber-1 and (i!=0 and i!=var_XNumber-1) and k==var_sliceNumber-1: # 情况13 下边界
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ##########四个边 end    ###########
                ############ 四个角 start ###############
                if i==0 and j==0 and k==var_sliceNumber-1: #########  情况14 原点O2点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==0 and k==var_sliceNumber-1: #########  情况15 A2点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==0 and j==var_YNumber-1 and k==var_sliceNumber-1: #########  情况16 B2点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==var_YNumber-1 and k==var_sliceNumber-1: #########  情况17 C2点
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k] + 2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ############ 四个角 end ###############
                if (i!=0 and i!=var_XNumber-1) and (j!=0 and j!=var_YNumber-1) and k==var_sliceNumber-1: #########  情况18 表面温度
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+2*d*(middle_temp_all[i][j][k-1]-middle_temp_all[i][j][k])
                ####### 二冷区出口处（z=L）网格节点的温度 end ################
                ####### 弯月面-----二冷区出口处之间 start ################
                ############## 四个边 start ############
                if i==0 and j==0 and (k!=0 and k!=var_sliceNumber-1):   #########  情况19 边O1-O2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==0 and j==var_YNumber-1 and (k!=0 and k!=var_sliceNumber-1):    #########  情况20 边B1-B2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==0 and (k!=0 and k!=var_sliceNumber-1):    #########  情况21 边A1-A2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and j==var_YNumber-1 and (k!=0 and k!=var_sliceNumber-1):    #########  情况22 边A1-A2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ############## 四个边 end ############
                ############## 四个面 start ############
                if i==0 and (j!=0 and j!=var_YNumber-1) and (k!=0 and k!=var_sliceNumber-1):  #########  情况23 面O1O2-B1B2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i+1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if (i!=0 and i!=var_XNumber-1) and j==0 and (k!=0 and k!=var_sliceNumber-1):  #########  情况24 面O1O2-A1A2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j+1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if (i!=0 and i!=var_XNumber-1) and j==var_YNumber-1 and (k!=0 and k!=var_sliceNumber-1):  #########  情况25 面C1C2-B1B2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+2*b*(middle_temp_all[i][j-1][k]-middle_temp_all[i][j][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltY*b*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                if i==var_XNumber-1 and (j!=0 and j!=var_YNumber-1) and (k!=0 and k!=var_sliceNumber-1):  #########  情况26 面C1C2-A1A2
                    h_int= h+EMISSIVITY*BOLTZMAN*(middle_temp_all[i][j][k]*middle_temp_all[i][j][k]+var_temperatureWater*var_temperatureWater)*(middle_temp_all[i][j][k]+var_temperatureWater)
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+2*a*(middle_temp_all[i-1][j][k]-middle_temp_all[i][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])-(2*h_int*deltX*a*(middle_temp_all[i][j][k]-var_temperatureWater))/Tconductivity
                ############## 四个面 end ############
                ####### 弯月面-----二冷区出口处之间 end ################
                ########### 不在任何边界上 start ###############
                if (i!=0 and i!=var_XNumber-1) and (j!=0 and j!=var_YNumber-1) and (k!=0 and k!=var_sliceNumber-1): #########  情况27  内部温度 
                    next_temp_all[i][j][k] = middle_temp_all[i][j][k]+a*(middle_temp_all[i+1][j][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i-1][j][k])+b*(middle_temp_all[i][j+1][k]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j-1][k])+d*(middle_temp_all[i][j][k+1]-2*middle_temp_all[i][j][k]+middle_temp_all[i][j][k-1])
                ########### 不在任何边界上 end ###############
    #### equation end #####
    return next_temp_all,sectionNumber
def select_h_in_different_section(k,var_ZNumber,var_Z,var_dis,var_h_initial,sectionNumber,var_VcastOriginal):
    delt_z=var_Z/var_ZNumber
    section_location_now=delt_z*k
    h=1.0
    if var_dis[0]<=section_location_now<=var_dis[1]:
        tTime = section_location_now/var_VcastOriginal
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
    return h,sectionNumber
def three_densonal_lifecycle_for_mpc(var_deltTime, var_XNumber, var_YNumber,var_X, var_Y,var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,var_latentHeatofSolidification,var_VcastOriginal,var_controlTime,var_sliceNumber,var_Z,var_dis,var_h_initial,var_castingTemp,delete_create_slice_number,slice_temp_all_initial,moveSliceNumber,stepTime_for_running):
    delt_Slice = var_Z/var_sliceNumber
    sectionNumber=[0]*(len(var_dis)-1) # 连铸二冷区各段切片数量
    next_temp_all = [([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    slice_temp_new =[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    middle_temp_all=[([([var_castingTemp]*var_sliceNumber)]*var_YNumber) for i in range(var_XNumber)]
    for k in range(var_sliceNumber):    
        for i in range(var_XNumber):
            for j in range(var_YNumber):
                slice_temp_new[i][j][k]=slice_temp_all_initial[i][j][k]
    mean_temp=[0]*(len(var_dis)-2) # 平均温度的初始化
    runningTimeinControltime = int(var_controlTime/var_deltTime)  #控制周期内的迭代次数
    sliceRunningTimeAll = int((delt_Slice/var_VcastOriginal)/var_deltTime) #切片间隔内运行的时间步长的次数
    centerPoint = int(var_YNumber/2)  # 铸坯宽度截面中心点的温度
    delt_delete_create_slice_number=moveSliceNumber-delete_create_slice_number
    #print(moveSliceNumber,delete_create_slice_number)
    for deltSliceNum in range(delt_delete_create_slice_number):
        for deltTime_Number in range(sliceRunningTimeAll):
            next_temp_all,sectionNumber=three_densonal_diff(var_deltTime,slice_temp_new,var_XNumber,var_YNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_latentHeatofSolidification,var_sliceNumber,var_Z,var_dis,var_h_initial,var_VcastOriginal)
        for k in range(var_sliceNumber):    
            for i in range(var_XNumber):
                for j in range(var_YNumber):
                    middle_temp_all[i][j][k]=next_temp_all[i][j][k]
        for k in range(var_sliceNumber):
            if k==0:
                for i in range(var_XNumber):
                    for j in range(var_YNumber):
                        slice_temp_new[i][j][0]=var_castingTemp
            else:
                for i in range(var_XNumber):
                    for j in range(var_YNumber):
                        slice_temp_new[i][j][k]=middle_temp_all[i][j][k-1]
    MoldNumber = sectionNumber[0]
    #for j in range(MoldNumber):
    #    print(next_temp_all[var_XNumber-1][centerPoint][j])
    sectionNumber_one = sectionNumber[0]+sectionNumber[1]
    sectionNumber_two = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]
    sectionNumber_three = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]
    sectionNumber_four = sectionNumber[0]+sectionNumber[1]+sectionNumber[2]+sectionNumber[3]+sectionNumber[4]
    print(sectionNumber)
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

    

