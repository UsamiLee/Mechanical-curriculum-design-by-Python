# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:38:16 2017

@author: UsamiHaru
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve

a = 0.310
L1 = 0.1424
L2 = 0.8627
L3 = 0.9379
L4 = 0.2958
L5 = 1.3317
b = 0.135
H = 0.250
K = 1.7
degF1 = 65 * math.pi / 180;
degF2 = 115 * math.pi / 180;
n2 = 65
w1 = 2 * np.pi *n2
plt.figure(1)
plt.figure(2)
plt.figure(3)
plt.figure(4)

VarDegree = np.linspace(0, 2 * math.pi, 360)


plt.figure(1)
plt.figure(figsize=(4,4))
coordinate_ax = L1 * np.cos(VarDegree)
coordinate_ay = L1 * np.sin(VarDegree)
plt.plot(coordinate_ax, coordinate_ay)
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.title('The displacement of Point A')
plt.show
   
def f(x):
    thet2, thet3, thet4 = x[0], x[1], x[2]
    a = x[3]
    return [L1 * np.cos(a) + L2 * np.cos(thet2) - L3 * np.cos(thet3) - L4 * np.cos(thet4),\
            L1 * np.sin(a) + L2 * np.sin(thet2) - L3 * np.sin(thet3) - L4 * np.sin(thet4),\
            np.sin(thet4) * L4 - b + L3 * np.sin(thet3),a - x[3]]

result = [] 
for i in range(len(VarDegree)):
    result.append( fsolve(f,[1.0,1.0,1.0,VarDegree[i]]) )    

Thet2 = []
Thet3 = []
Thet4 = []    
for  i in range(len(VarDegree)):
    Thet2.append(result[i][0])
    Thet3.append(result[i][1])
    Thet4.append(result[i][2])
Thet2 = np.array(Thet2)
Thet3 = np.array(Thet3)
Thet4 = np.array(Thet4)
    
plt.figure(2)
plt.figure(figsize=(10,4))
coordinate_bx = VarDegree
angularvelocityw3 = (L1 * np.sin(VarDegree - Thet2) * w1) / (L3 - (L3 - L4) * np.sin(Thet3 - Thet2)) 
c = np.polyfit(coordinate_bx, angularvelocityw3,3)
yy = np.polyval(c,coordinate_bx)
f_liner =np.polyval(c,coordinate_bx)
plt.plot(coordinate_bx, f_liner)
plt.xlabel('The original arc')
plt.ylabel('rad/s')
plt.title('The angularvelocity of Shaft DE')
plt.show

plt.figure(3)
plt.figure(figsize = (10,4))
coordinate_cx = VarDegree
angularvelocityw4 = - ((angularvelocityw3 * L3 * np.cos(Thet3)) / (L4 * np.cos(Thet4)))
sildervelocity = - angularvelocityw3 * L3 * np.sin(Thet3) - L4 * np.sin(Thet4) * angularvelocityw4
c_1 = np.polyfit(coordinate_cx, sildervelocity,3)
yy_1 = np.polyval(c_1,coordinate_cx)
f_liner_1 = np.polyval(c_1,coordinate_cx)                                                 
plt.plot(coordinate_cx,f_liner_1)
plt.xlabel('The original arc')
plt.ylabel('m/s')
plt.title('The velocity of Slider')
plt.show

w3 = angularvelocityw3
w4 = angularvelocityw4
Lce = L3 - L4
w2 = (w3 * (L3 - Lce) * np.sin(Thet3) - w1 * L1 * np.sin(VarDegree)) / (L2 * np.sin(Thet2))
D = (L3 - Lce) * w3 * w3 * np.cos(Thet3) - L1 * (w1 * w1 * np.cos(VarDegree)) - L2 * w2 * w2 * np.cos(Thet2)
E = -(L3 - Lce) * w3 * w3 * np.sin(Thet3) - L1 *( -w1 * w1 * np.sin(VarDegree) + L2 * w2 * w2 * np.sin(Thet2))
e3 = (D * np.cos(Thet2) - E * np.sin(Thet2)) / (L3 - Lce) * np.sin(Thet3 - Thet2)
e2 = (D + (L3 - Lce) * e3 * np.sin(Thet3)) / (L2 * np.sin(Thet2))
e4 = (L4 * w4 * w4 * np.cos(Thet4) - L3 * (e3 * np.cos(Thet3) - w3 * w3 * np.sin(Thet3))) / (L4 * np.cos(Thet4))
plt.figure(4)
plt.figure(figsize = (10,4))
af = -L3 * (e3 * np.sin(Thet3) + w3 * w3 * np.cos(Thet3)) - L4 * (e4 * np.sin(Thet4) + w4 * w4 * np.cos(Thet4))
c_2 = np.polyfit(VarDegree, af, 3)
yyy = np.polyval(c_2,VarDegree)
f_liner_2 = np.polyval(c_2, VarDegree)
plt.plot(VarDegree,f_liner_2)
plt.xlabel('The original arc')
plt.ylabel('m/s^2')
plt.title('The acceleration of Slider')
plt.show