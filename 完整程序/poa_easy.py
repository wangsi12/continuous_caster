# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 19:32:22 2019

@author: ws
"""
#coding=utf-8

import random
from GetResult import ObjectiveFunction
#随机初始果蝇群体位置
def smell_sum(s):
    s,result=ObjectiveFunction(s)
    return result

#print(x_init,y_init)
maxgen=10 #迭代次数
sizepop=20 #种群规模
variable_size=4#决策变量个数
x_init=[1.0*random.random() for i in range(variable_size)]
y_init=[1.0*random.random() for i in range(variable_size)]
x=[[0 for i in range(variable_size)] for j in range(sizepop)]
y=[[0 for i in range(variable_size)] for j in range(sizepop)]
D=[[0 for i in range(variable_size)] for j in range(sizepop)]
S=[[0 for i in range(variable_size)] for j in range(sizepop)]
Smell=[0 for i in range(sizepop)]
#每代最优Smell值纪录到yy数组中，并记录最优迭代坐标
yy=[0 for i in range(maxgen)]
X_best=[0 for i in range(maxgen)]
Y_best=[0 for i in range(maxgen)] 
#果蝇寻优开始，利用嗅觉寻找食物
for i in range(sizepop):
    for j in range(variable_size):
        x[i][j]=x_init[j]+0.2*random.random()-0.1
        y[i][j]=y_init[j]+0.2*random.random()-0.1
        D[i][j]=(x[i][j]**2+y[i][j]**2)**(0.5)
        S[i][j]=100/D[i][j]
    Smell[i]=smell_sum(S[i])#记录每代所有种群的适应值
#找出初始最是适应值
bestSmell=min(Smell)
bestindex=Smell.index(bestSmell)
#保留最佳味道浓度值与x、y坐标，此时果蝇群体利用视觉往该位置飞去。
x_init=x[bestindex]
y_init=y[bestindex]
Smellbest=bestSmell
#果蝇迭代寻优开始
for g in range(maxgen):
#附与果蝇个体利用嗅觉搜寻食物之随机方向与距离
    for i in range(sizepop):
        for j in range(variable_size):
            x[i][j]=x_init[j]+2*random.random()-1
            y[i][j]=y_init[j]+2*random.random()-1
            D[i][j]=(x[i][j]**2+y[i][j]**2)**(0.5)
            S[i][j]=100/D[i][j]
        Smell[i]=smell_sum(S[i])#记录每代所有种群的适应值

    #找出此果蝇群体的中味道浓度最高的果蝇(求极大值)
    bestSmell=min(Smell)
    bestindex=Smell.index(bestSmell)
    print("bestSmell",bestSmell)
    print("S",S)
    #判断味道浓度是否优于前一迭代味道浓度，
    # 若是则保留最佳味道浓度值与x、y坐标，此时果蝇群体利用视觉往该位置飞去。
    if bestSmell<Smellbest:
        x_init=x[bestindex];y_init=y[bestindex]
        Smellbest=bestSmell

#每代最优Smell值纪录到yy数组中，并记录最优迭代坐标
    yy[g]=Smellbest
    
    X_best[g]=x_init;Y_best[g]=y_init
    #print(S)
print("")

print ("经过果蝇算法优化后的最优结果：")
print(yy)
print (Smellbest)



