from one_dimensional_mpc_solver import one_dimensional_mpc_solver
var_XNumber = 30
var_ZNumber = 1500 # 拉坯方向网格划分的数量
var_Z = 31.332     # 拉坯长度
deltZ = var_Z/var_ZNumber # 拉皮方向的空间间隔
middle_temp = [([1530]*var_ZNumber) for i in range(var_XNumber)]
var_X = 0.32
var_Y=1.25 # 铸坯宽度
var_temperatureWater = 30
var_rouS=7800 
var_rouL=7200 
var_specificHeatS=660 
var_specificHeatL=830 
var_TconductivityS=31 
var_TconductivityL=35 
var_liqTemp=1514 
var_SodTemp=1453 
var_m=1.34 
var_latentHeatofSolidification = 268000
var_VcastOriginal=0.6/60
var_controlTime=5
var_dis=[0.0,0.9,3.86,5.679,7.574,11.363,17.457,23.402,29.022,31.332] # 连铸二冷区各段距弯月面的距离 单位m
var_runningTime=50
var_castingTemp=1530 # 铸坯的浇筑温度
var_aimTemp=[1000,985,942,941,915,860,812,800] # 测量的温度值
var_SCZ_num=8 # 连铸二冷区冷却区段的数量
Var_runningSteadyTime = 10 #稳态运行时间
var_VcastNew=1.0/60 # 拉速
var_runningTime = 100 # 运行时间
var_controlTime = 5 # 控制时间
######### 需要输入的参数 end   #################

#############  parameters for optimation start ################
var_deltTime=0.1 # 差分计算时间间隔
var_iter_max=3; # 优化迭代的最大次数
var_eps=200 # 优化的停止标准
var_water_initial = [300,260,240,200,180,160,100,80] # 二冷水量初值
var_alfa=[545.4,381.4,246.4,182.5,158.3,147.9,147.8,146.9] # 待定参数的初值
#############  parameters for optimation end ################
h_end_all=[]
h_end_all=one_dimensional_mpc_solver(var_controlTime,var_runningTime,var_ZNumber,var_SCZ_num,var_XNumber,var_X,var_Y,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_aimTemp,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_water_initial,var_latentHeatofSolidification,var_ZNumber,var_Z,var_VcastNew,var_alfa)
