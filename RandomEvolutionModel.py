#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 19:07:01 2021

@author: wyattpetryshen
"""

# Code templates for Ornstein-Uhlenbeck process and Brownian motion are from IPython Interactive Computing and Visualization Cookbook, Second Edition (2018), by Cyrille Rossant.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import stats
import time as timetime
import random
import itertools


#Calculate angle between vectors
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def vectorTransform(end_Point,old_origin):
    """ Returns vector translated to origin."""
    newP = np.subtract(end_Point,old_origin)
    return newP

def angle_between(v1, v2):
    """ Returns the angle in degrees between vectors 'v1' and 'v2'."""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
#Equation for Ornstein-Uhlenbeck process

def binVals(array,step):
    """ Bins array and calculates mean for specified bin size."""
    start = 0
    step = step
    stop = step
    iterations = int(len(array)/step)
    meanvals = []
    for i in np.arange(0,iterations):
        tempbin = array[start:stop]
        meanbin = np.mean(tempbin)
        start = start + step
        stop = stop + step
        meanvals.append(meanbin)
    return(meanvals)

#Code for figure 1 in supplementary information
#Change the parameters accordingly

#Sine wave
sample_rate = 1000
time = np.arange(0, 10, 1/sample_rate)
frequency = 0.1
amplitude = 4
theta = 0
sinewave = amplitude * np.sin(2 * np.pi * frequency * time + theta)

##Model parameters
sigma = 1  #standard deviation
mu = 0  #mean
tau = 0.05  #time constant

##simulation parameters
dt = 0.0001  #Time step
T = 1  #Total time
n = int(T/dt)  #Number of steps
t = np.linspace(0., T, n)  #Vector of times

##Calculated randomized variables
sigma_bis = sigma * np.sqrt(2. / tau)
sqrtdt = np.sqrt(dt)

#Plot of Sine wave
plt.plot(time,sinewave)
plt.title(r'SineWave with amplitude = {}, frequency = {}'.format(amplitude,frequency))
plt.axis([0, 10, -4, 4])

#Random Drift
for iters in range(100):
    ##Store results
    x = np.zeros(n)
    #Euler-Maruyama method
    for i in range(n - 1):
        x[i + 1] = x[i] + dt * (-(x[i] - sinewave[i]) / tau) + sigma_bis * sqrtdt * np.random.randn()
    array = x
    plt.plot(time,array,linewidth=0.1)
    plt.title(r'OH Drift with amplitude = {}, frequency = {}'.format(amplitude,frequency))

#Time-averaged drift
for iters in range(100):
    ##Store results
    x = np.zeros(n)
    #Euler-Maruyama method
    for i in range(n - 1):
        x[i + 1] = x[i] + dt * (-(x[i] - sinewave[i]) / tau) + sigma_bis * sqrtdt * np.random.randn()
    array = x
    meanarray = binVals(array,int(sample_rate))
    plt.plot(time[int(sample_rate/2):-1:int(sample_rate)],meanarray,linewidth=0.1)
    plt.scatter(time[int(sample_rate/2):-1:int(sample_rate)],meanarray,linewidth=0.1)
    plt.title(r'OH Drift time-averaged with amplitude = {}, frequency = {}'.format(amplitude,frequency))
    plt.axis([0, 10, -6, 6])
    #plt.plot(time,x)


#Iterate OH means and calculate the angle between vectors
start_time = timetime.time()
angle_list = []
for iters in range(100):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    for i in range(n - 1):
        x1[i + 1] = x1[i] + dt * (-(x1[i] - sinewave[i]) / tau) + sigma_bis * sqrtdt * np.random.randn()
        x2[i + 1] = x2[i] + dt * (-(x2[i] - sinewave[i]) / tau) + sigma_bis * sqrtdt * np.random.randn()
    meanarray1, meanarray2 = binVals(x1,int(sample_rate)),binVals(x2,int(sample_rate))
    for j in np.arange(1,len(meanarray1)):
        if j != len(meanarray1)-1:
            Idx_O = j
            Idx_E = j+1
            v1 = vectorTransform((time[int(sample_rate/2):-1:int(sample_rate)][Idx_O],meanarray1[Idx_O]),(time[int(sample_rate/2):-1:int(sample_rate)][Idx_E],meanarray1[Idx_E]))
            v2 = vectorTransform((time[int(sample_rate/2):-1:int(sample_rate)][Idx_O],meanarray2[Idx_O]),(time[int(sample_rate/2):-1:int(sample_rate)][Idx_E],meanarray2[Idx_E]))
            vector_angle = angle_between(v1,v2)
            angle_list.append(vector_angle)
        else:
            pass
plt.hist(angle_list, bins = np.arange(0,180,5))
plt.xlabel('Angle')
plt.ylabel('Probability')
plt.title(r'Histogram of OH Trait Drift a=4, frequency=0.1')
print("--- %s seconds ---" % (timetime.time() - start_time))

###Brownian motion
#simulation parameters
n = 100000 #time steps

#Two one dimensional cases that can be combined into two dimensional case
x = np.cumsum(np.random.randn(n))
y = np.cumsum(np.random.randn(n))

xP = np.cumsum(np.random.randn(n))
yP = np.cumsum(np.random.randn(n))

# We add 10 intermediary points between two
# successive points. We interpolate x and y.
k = 50

x2 = np.interp(np.arange(n * k), np.arange(n) * k, x)
y2 = np.interp(np.arange(n * k), np.arange(n) * k, y)
xP2 = np.interp(np.arange(n * k), np.arange(n) * k, xP)
yP2 = np.interp(np.arange(n * k), np.arange(n) * k, yP)

# Now, we draw our points with a gradient of colors.
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
plt.scatter(x2, y2, c=range(n * k), linewidths=0,
           marker='o', s=3, cmap=plt.cm.jet,)
plt.axis('equal')

plt.scatter(xP2, yP2, c = range(n*k), linewidths=0,
           marker='o', s=3, cmap=plt.cm.jet,)

plt.plot(xP, yP)
plt.plot(x2, y2)


#Brownian Time-averaged drift for single lineages
for iters in range(1000):
    ##Store results
    n = 1000 #time steps
    #Two one dimensional cases that can be combined into two dimensional case
    x = np.cumsum(np.random.randn(n))
    y = np.cumsum(np.random.randn(n))
    #Find mean values
    meanx = binVals(x,int(100))
    meany = binVals(y,int(100))
    #plot
    plt.plot(meanx,meany,linewidth=0.5)
    plt.scatter(meanx,meany,linewidth=0.1, s = 4)
    plt.title(r'Brownian Drift time-averaged with equal rates')
    #plt.axis([0, 10, -6, 6])


#Iterate BM means and calculate the angle between vectors
start_time = timetime.time()
BW_angle_list = []
for iters in range(10000):
    runs = 1000 #time steps
    Tavg = 100 #Average years; for random number use random.randrange(i,j)
    rate_One = 1 #rate change
    rate_Two = 10 #rate change
    x, y = np.cumsum(np.random.randn(runs * rate_One)), np.cumsum(np.random.randn(runs * rate_One))
    x2, y2 = np.cumsum(np.random.randn(runs * rate_Two)), np.cumsum(np.random.randn(runs * rate_Two))
    meanx,meany = binVals(x,int(Tavg)*rate_One),binVals(y,int(Tavg)*rate_One)
    meanx2,meany2 = binVals(x2,int(Tavg)*rate_Two),binVals(y2,int(Tavg)*rate_Two)
    for j in np.arange(1,len(meanx)):
        if j != len(meanx)-1:
            Idx_O = j
            Idx_E = j+1
            v1 = vectorTransform((meanx[Idx_O],meany[Idx_O]),(meanx[Idx_E],meany[Idx_E]))
            v2 = vectorTransform((meanx2[Idx_O],meany2[Idx_O]),(meanx2[Idx_E],meany2[Idx_E]))
            vector_angle = angle_between(v1,v2)
            BW_angle_list.append(vector_angle)
        else:
            pass
plt.hist(BW_angle_list, bins = np.arange(0,180,1))
plt.xlabel('Angle')
plt.ylabel('Probability')
plt.title(r'Histogram of BW Parallelism')
print("--- %s seconds ---" % (timetime.time() - start_time))
