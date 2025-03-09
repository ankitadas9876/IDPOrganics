#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import csv
import seaborn as sns
from numpy import loadtxt
import pandas as pd
import sys
import os
from subprocess import call
import math as mt
from numpy import trapz
import time
import glob
FlynnIN = [2.595*10**-11, 2.725*10**-10, 9.343*10**-10, 8.954*10**-10, 5.840*10**-10, 2.336*10**-10, 9.084*10**-11, 7.786*10**-12, 2.595*10**-12, 1.168*10**-12, 3.893*10**-13, 7.786*10**-14]
FlynnIN = (np.array(FlynnIN,dtype=np.float32))
print(len(FlynnIN))
Mparticle = [3.94569*10**-7, 4.01944*10**-8, 3.88242*10**-9, 3.81704*10**-10, 4.30989*10**-11, 3.88242*10**-12, 3.81704*10**-13, 4.77129*10**-14, 4.18879*10**-15, 2.68083*10**-16, 4.77129*10**-17, 4.18879*10**-18 ]
IDPdiam1 = np.array([910.0, 425.0, 195.0, 90.0, 43.5, 19.5, 9, 4.5, 2.0, 0.8, 0.45, 0.2],dtype=np.float32)
IDPdiam = IDPdiam1/(1000000)
Surfload = [None]*len(IDPdiam)
IDPunalt = np.array([0.01, 0.07, 0.24, 0.46, 0.75, 0.88, 0.95, 1, 1, 1, 1, 1])
IDPunmelt = np.array([0.5, 0.7, 0.9, 0.98, 1, 1, 1, 1, 1, 1, 1, 1])
path= "C://Users/arind/Downloads/"
sys.path.append(path)
dust = pd.DataFrame(np.loadtxt('normalizeddust.txt'))
dust = np.array(dust)
dust_new = np.zeros((47,63))
for i in range(47):
    for j in range(63):
        dust_new[i,j] = dust[i+1,j+1]
print(len(dust_new))
fig = plt.figure()
fig.set_size_inches(11.5, 9.5)
im = plt.imshow((10**6*dust_new/(870)), cmap=cm.RdBu, extent=(-177, 177, -88.1 , 88.1))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Normalized Global Dust Deposition From MCD data')
plt.gca().invert_yaxis()


plt.colorbar(im, label ='Dust deposited in $u$m')


# In[3]:


dfc = pd.DataFrame(np.loadtxt('uvcdata.txt'), columns=['long', 'lat', 'energyc'])
dfb = pd.DataFrame(np.loadtxt('uvbdata.txt'), columns=['long', 'lat', 'energyb'])
dfa = pd.DataFrame(np.loadtxt('uvadata.txt'), columns=['long', 'lat', 'energya'])
dfd = pd.DataFrame(np.loadtxt('tempdata.txt'), columns=['long', 'lat', 'temp'])
lat = np.array(dfc['lat'])
energy_c = np.array(dfc['energyc'])
energy_b = np.array(dfb['energyb'])
energy_a = np.array(dfa['energya'])
temp=np.array(dfd['temp'])
finalcarbon = np.ones((47,63),dtype=np.float32)
percarbon = np.ones((47,63),dtype=np.float32)
for j in range(1,48):
    e_new_c = [None]*36
    e_totalc = 0
    e_new_b = [None]*36
    e_totalb = 0
    e_new_a = [None]*36
    e_totala = 0
    temp_new = [None]*36
    k=0
    for i in range(len(energy_c)):
        if(lat[i]==(dust[j,0]//10)*10):
            x = [lat[i],lat[i+1]]
            y1 = [energy_c[i],energy_c[i+1]]
            y2 = [energy_b[i],energy_b[i+1]]
            y3 = [energy_a[i],energy_a[i+1]]
            y4 = [temp[i],temp[i+1]]
            e_new_c[k] = np.interp(dust[j,0],x,y1)
            e_totalc = e_totalc + e_new_c[k]
            e_new_b[k] = np.interp(dust[j,0],x,y2)
            e_totalb = e_totalb + e_new_b[k]
            e_new_a[k] = np.interp(dust[j,0],x,y3)
            e_totala = e_totala + e_new_a[k]
            temp_new[k] = np.interp(dust[j,0],x,y4)
            k = k+1
    for i in range(1,64):
        #print(dust[j,0])
        #print(dust[0,i])
        years = np.array(np.linspace(1,10,num=10))
        Ls = np.array(list(range(0,360,10)))
        #energy_a = e_new_a
       # energy_b = e_new_b
       # energy_c = e_new_c
        energy_time_c = np.array(np.zeros((len(years),len(Ls))),dtype=np.float32)
        energy_time_b = np.array(np.zeros((len(years),len(Ls))),dtype=np.float32)
        energy_time_a = np.array(np.zeros((len(years),len(Ls))),dtype=np.float32)
        IDPlife = np.array(np.zeros((len(years),len(IDPdiam))),dtype=np.float32)
        etotal_Ls = np.array(np.zeros((len(years),len(Ls))))
        ec = np.array(np.zeros(len(years)),dtype=np.float32)
        eb = np.array(np.zeros(len(years)),dtype=np.float32)
        ea = np.array(np.zeros(len(years)),dtype=np.float32)
        ec[0] = sum(e_new_c)
        eb[0] = sum(e_new_b)
        ea[0] = sum(e_new_a)
        energy_time_c[0,:] = e_new_c
        energy_time_b[0,:] = e_new_b
        energy_time_a[0,:] = e_new_a
        rate = dust[j,i]/870 #density of JSC-1 Mars 0.87 kg/cm3 (Allen et al., 1998)
        
        for l in range(len(years)):
            for m in range(len(Ls)):
                
                energy_time_c[l,m] = energy_c[m]/(2.71**(l*30.49*1000*rate))
                ec[l] = ec[l] + energy_time_c[l,m]
                energy_time_b[l,m] = energy_b[m]/(2.71**(l*29.54*1000*rate))
                eb[l] = eb[l] + energy_time_b[l,m]
                energy_time_a[l,m] = energy_a[m]/(2.71**(l*29.35*1000*rate))
                ea[l] = ea[l] + energy_time_a[l,m]
                
        QE = np.array([None]*len(Ls),dtype=np.float32)
        FQEdaily = np.array(np.zeros((len(years),len(Ls)),dtype=np.float32))
        totalfqe = np.array([None]*len(years),dtype=np.float32)
        totaluv = np.array(np.zeros(len(years)),dtype=np.float32)
        IDPlife = np.array(np.zeros((len(years),len(IDPdiam)),dtype=np.float32))
        for l in range(len(Ls)):
            QE[l] = (temp_new[l]*((2.76*10**-15))-(1.54*10**-13)) 
            
        for l in range(len(years)):
            for m in range(len(Ls)):
                totaluv[l] = (energy_time_c[l,m]+energy_time_b[l,m]+energy_time_a[l,m])
                FQEdaily[l,m] = QE[l]*(energy_time_c[l,m]+energy_time_b[l,m]+energy_time_a[l,m])
       
        for l in range(len(years)):
            totalfqe[l] = 0.0
            for m in range(len(Ls)):
                totalfqe[l] = totalfqe[l]+FQEdaily[l,m]
    
        totalfqe = (totalfqe/(88775.0*36.0))*(59354500.0) #88775 = seconds in one martian year, 36= number of solar longitude points,59354500 = FUV average
        for l in range(len(years)):
            for m in range(len(IDPdiam)):
                IDPlife[l,m] = (2*0.1*IDPunalt[m]*1000*IDPunmelt[m]*IDPdiam[m])/(3*totalfqe[l]*0.012) 
           
        carbon = np.array(np.zeros((len(years),len(IDPdiam)),dtype=np.float32))
        carbon[0,:] = (np.array(FlynnIN))
        sumcarbon = np.array([None]*len(years),dtype=np.float32)
        layermass = np.zeros(len(years),dtype=np.float32)
        ppb = np.zeros(len(years),dtype=np.float32)
        for l in range(len(years)):
            layermass[l] = rate*870 #density of JSC-1 Mars 0.87 kg/cm3 (Allen et al., 1998)  
        for l in range(len(years)-1):
            for m in range(len(IDPdiam)):
                carbon[l+1,m] = (carbon[l,m] - carbon[l,m]/(IDPlife[l,m]))
        for l in range(len(years)):
            sumcarbon[l] = sum(carbon[l,:])
            ppb[l] = (sumcarbon[l]/layermass[l])*10**9
        if(ppb[l]<10**4 and ppb[l]>100):
            finalcarbon[j-1,i-1]=ppb[l]
            percarbon[j-1,i-1] = (ppb[0]-ppb[l])/ppb[l]
        else:
            finalcarbon[j-1,i-1] = None
            percarbon[j-1,i-1] = None
        


# In[4]:


import cmcrameri.cm as cmc
import scienceplots 
plt.style.use(['science','notebook','grid'])
'''fig = plt.figure()
fig.set_size_inches(11.5, 9.5)
im = plt.imshow((np.log10(finalcarbon)), cmap=cm.RdBu, extent=(-177, 177, -88.1 , 88.1))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Organics Preserved Globally using Normalized dust deposition')
plt.gca().invert_yaxis()
plt.colorbar(im, label ='Organics preserved in ppb (log scale)')'''

plt.rc("font", family="serif", size=18.)
fig = plt.figure()
fig.set_size_inches(11.5, 9.5)
im = plt.imshow(((100-percarbon*100)), cmap=cmc.batlow, extent=(-177, 177, -88.1 , 88.1))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('% Organics Preserved Globally using Normalized dust deposition')
plt.gca().invert_yaxis()
plt.colorbar(im, label ='percentage preserved after 10 Martian years')
plt.savefig('IDPOrganicMarsMap_Das')


# In[ ]:





# In[ ]:




