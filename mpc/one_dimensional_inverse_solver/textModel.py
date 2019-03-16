from one_dimensional_inverse_solver import one_dimensional_inverse_solver
#############  parameters for continuous casting start ################
var_XNumber=30 # 铸坯厚度方向的网格点数
var_X=0.32 # 铸坯厚度
var_Y=1.25 # 铸坯宽度
var_temperatureWater=30 # 冷却水的温度
var_rouS=7800 # 钢液凝固以后的密度
var_rouL=7200 # 液态钢水的密度
var_specificHeatS=660 # 钢液凝固以后的比热
var_specificHeatL=830 # 液态钢水的比热
var_TconductivityS=31 # 钢液凝固以后的导热系数
var_TconductivityL=35 # 液态钢水的导热系数
var_liqTemp=1514 # 钢的液相线温度
var_SodTemp=1453 # 钢的固相线温度
var_m=1.34 # 常数
var_castingTemp=1544 # 铸坯的浇筑温度
var_dis=[0.0,0.9,3.86,5.679,7.574,11.363,17.457,23.402,29.022,31.332] # 连铸二冷区各段距弯月面的距离 单位m
var_NumberofSCZ = 8
var_measuredTemp=[1010,985,942,941,915,860,812] # 测量的温度值
var_disMeasure=[3.9,5.8,11.3,13.5,17.9,23.5,31] # 测温点距弯月面的距离
var_numberofMeasuredPoint=7 # 测量点的个数
var_SCZ_num=8 # 连铸二冷区冷却区段的数量
var_VcastOriginal=0.6/60 # 铸坯的拉坯速度 m/min
latentHeatofSolidification = 26800 #凝固潜热 
#############  parameters for continuous casting end ################

#############  parameters for optimation start ################
var_deltTime=0.1 # 差分计算时间间隔
var_iter_max=10; # 优化迭代的最大次数
var_eps=200 # 优化的停止标准
var_alfa_initial=[4,5,4.5,3.6,3,2.6,4,3.8] # 待定参数的初值
var_WaterFlowRate = [300,260,240,200,180,160,100,80]
#############  parameters for optimation end ################
alfa_get=[]
alfa_get=one_dimensional_inverse_solver(var_XNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_measuredTemp,var_disMeasure,var_numberofMeasuredPoint,var_SCZ_num,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_alfa_initial,latentHeatofSolidification,var_WaterFlowRate,var_NumberofSCZ)
