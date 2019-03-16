# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 20:00:18 2019

@author: ws
"""

import math
import random
from sklearn.datasets import load_digits
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
#from sklearn.neural_network import MLPClassifier
import os
import sys
import time
 
import numpy
 
import theano
import theano.tensor as T

import xlrd
import xlwt
import matplotlib.pyplot as plt

worksheet = xlrd.open_workbook('data800.xlsx')
sheet_names= worksheet.sheet_names() 

sheet1 = worksheet.sheet_by_name('Sheet1')
sheet2 = worksheet.sheet_by_index(0)
v = sheet2.col_values(1)
del v[0]
v2 = sheet2.col_values(2)
del v2[0]
T = sheet2.col_values(3)
del T[0]
Q1 = sheet2.col_values(4)
del Q1[0]
Q2 = sheet2.col_values(5)
del Q2[0]
Q3 = sheet2.col_values(6)
del Q3[0]

shuru=[]
for i in range(len(v)):
    row=[]
    row.append(v[i])
    row.append(v2[i])
    row.append(T[i])
    
    shuru.append(numpy.array(row))
#print(shuru)
q1=[]
for i in range(len(Q1)):
    q1.append(int(Q1[i]))
q2=[]
for i in range(len(Q2)):
    q2.append(int(Q2[i]))
q3=[]
for i in range(len(Q3)):
    q3.append(int(Q3[i]))

shucu=[]
for i in range(len(v)):
    row=[]
    row.append(q1[i])
    row.append(q2[i])
    row.append(q3[i])
    
    shucu.append(numpy.array(row))
#Q1    
X = numpy.array(shuru)
y = numpy.array(q1)
#digits = load_digits()
#X = digits.data
#y = digits.target

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=3)
random.seed(0)
print("y_train",y_train)
print(y_test)

Y_train=[]
for i in y_train:
    s=[i]
    Y_train.append(s)
    
import numpy as np
 
#定义双曲函数和他们的导数
def tanh(x):
    return np.tanh(x)
 
def tanh_deriv(x):
    return 1.0 - np.tanh(x)**2
 
def logistic(x):
    return 1/(1 + np.exp(-x))
 
def logistic_derivative(x):
    return logistic(x)*(1-logistic(x))
def relu(x): 
    
    return np.maximum(x, 0.0)
def relu_derivative(x):
    return 1
 
#定义NeuralNetwork 神经网络算法
class NeuralNetwork:
    #初始化，layes表示的是一个list，eg[10,10,3]表示第一层10个神经元，第二层10个神经元，第三层3个神经元
    def __init__(self, layers, activation='relu'):
        if activation == 'logistic':
            self.activation = logistic
            self.activation_deriv = logistic_derivative
        elif activation == 'tanh':
            self.activation = tanh
            self.activation_deriv = tanh_deriv
        elif activation == 'relu':
            self.activation = relu
            self.activation_deriv = relu_derivative
        print(len(layers) - 1)
        self.weights = []
        #循环从1开始，相当于以第二层为基准，进行权重的初始化
        for i in range(1, len(layers) - 1):
            #对当前神经节点的前驱赋值
            self.weights.append((2*np.random.random((layers[i - 1] + 1, layers[i] + 1))-1)*0.25)
            #对当前神经节点的后继赋值
            #numpy.random.rand(d0,d1,…,dn)dn表示每个维度,这里表示layers[i - 1] + 1个layers[i] + 1维数
            #self.weights[0]是一个完整的array
            self.weights.append((2*np.random.random((layers[i] + 1, layers[i + 1]))-1)*0.25)
        #print(self.weights[0])
    #训练函数   ，X矩阵，每行是一个实例 ，y是每个实例对应的结果，learning_rate 学习率， 
    # epochs，表示抽样的方法对神经网络进行更新的最大次数
    def fit(self, X, y, learning_rate=0.00001, epochs=10000):
        X = np.atleast_2d(X) #确定X至少是二维的数据
           
        temp = np.ones([X.shape[0], X.shape[1]+1]) #初始化矩阵
        temp[:, 0:-1] = X  # adding the bias unit to the input layer
        X = temp
        print(len(X))
        y = np.array(y) #把list转换成array的形式
     
        for k in range(50*len(X)):
            #随机选取一行，对神经网络进行更新
            #i = np.random.randint(X.shape[0]) 
            i=k%50
            a = [X[i]]
            #print(a)
            #完成所有正向的更新
            
            for l in range(len(self.weights)):
                a.append(self.activation(np.dot(a[l], self.weights[l])))
                #print(self.activation(10))
            #
            error = y[i] - a[-1]
            #print( a[-1])
            deltas = [error * self.activation_deriv(a[-1])]
 
            #开始反向计算误差，更新权重
            for l in range(len(a) - 2, 0, -1): # we need to begin at the second to last layer
                deltas.append(deltas[-1].dot(self.weights[l].T)*self.activation_deriv(a[l]))
            deltas.reverse()
            for i in range(len(self.weights)):
                layer = np.atleast_2d(a[i])
                delta = np.atleast_2d(deltas[i])
                self.weights[i] += learning_rate * layer.T.dot(delta)
    
    #预测函数            
    def predict(self, x):
        x = np.array(x)
        #print(x)
        temp = np.ones(x.shape[0]+1)
        #print(temp)
        temp[0:-1] = x
        a = temp
        for l in range(0, len(self.weights)):
            #print(np.dot(a, self.weights[l]))
            a = self.activation(np.dot(a, self.weights[l]))
           
        return a    
import numpy as np

nn = NeuralNetwork([3,10,1], 'relu')
#X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
X=X_train
#y = np.array([0, 1, 1, 0])
y=y_train
#print(X)
#print(y)
print(y_test)
y_predict=[]
nn.fit(X, y)
for i in X_test:
    y_predict.append(int(nn.predict(i)))
print(y_predict)
y_ad=0
for i in range(len(X_test)):
    y_ad+=np.sqrt((y_test[i]-y_predict[i])**2)
    #print(y_ad)
print(y_ad/len(y_test))
print("Q1预测结果图")
x = range(len(y_predict))
y=y_test
y2=y_predict
plt.figure(figsize=(100,4))
plt.plot(x,y)
plt.plot(x,y2,color='red')
#plt.plot(x,y3,color='black',linestyle='--')
plt.show()
x = range(len(y_predict))


y=y_test
y2=y_predict
plt.figure(figsize=(100,4))
plt.scatter(x,y)
plt.scatter(x,y2,color='red')
#plt.plot(x,y3,color='black',linestyle='--')
plt.show()