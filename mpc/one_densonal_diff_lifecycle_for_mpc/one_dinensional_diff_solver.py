def one_densonal_diff(h, var_deltTime, middle_temp, var_XNumber, var_X, var_temperatureWater, var_rouS, var_rouL, var_specificHeatS, var_specificHeatL, var_TconductivityS, var_TconductivityL, var_liqTemp, var_SodTemp, var_m,latentHeatofSolidification):
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
                var_liqTemp-var_SodTemp)+var_specificHeatL+latentHeatofSolidification/(var_liqTemp-var_SodTemp)
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
