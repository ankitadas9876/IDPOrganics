#!/usr/bin/env python
# coding: utf-8

# In[17]:


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
import scienceplots 
import math
from math import e
plt.style.use(['science','notebook','grid'])
path = "C://Users/arind/Downloads"
sys.path.append("C://Users/arind/Downloads")
plt.rcParams.update({'font.size': 10})

'''Script to calculate remaining Organic at Gale 
dust deposition value: 3.5 microns
Latitude: -4.5
'''

#FlynnIN: Incoming total IDP organic infall in kg per m² per Martian Year
FlynnIN = [2.595*10**-11, 2.725*10**-10, 9.343*10**-10, 8.954*10**-10, 5.840*10**-10, 2.336*10**-10, 9.084*10**-11, 7.786*10**-12, 2.595*10**-12, 1.168*10**-12, 3.893*10**-13, 7.786*10**-14]
Mparticle = [3.94569*10**-7, 4.01944*10**-8, 3.88242*10**-9, 3.81704*10**-10, 4.30989*10**-11, 3.88242*10**-12, 3.81704*10**-13, 4.77129*10**-14, 4.18879*10**-15, 2.68083*10**-16, 4.77129*10**-17, 4.18879*10**-18 ]
IDPdiam1 = np.array([910, 425, 195, 90, 43.5, 19.5, 9, 4.5, 2, 0.8, 0.45, 0.2]) #IDP particle diameters in micron
IDPdiam = IDPdiam1/(10**6) #IDP diameter in meters
Surfload = [None]*len(IDPdiam)
IDPunalt = np.array([0.01, 0.07, 0.24, 0.46, 0.75, 0.88, 0.95, 1, 1, 1, 1, 1]) #IDP fraction un-altered for each particle size
IDPunmelt = np.array([0.5, 0.7, 0.9, 0.98, 1, 1, 1, 1, 1, 1, 1, 1]) #IDP fraction un-melted for each particle size


# In[18]:


dfc = pd.DataFrame(np.loadtxt('uvcdata.txt'), columns=['long', 'lat', 'energyc'])
dfb = pd.DataFrame(np.loadtxt('uvbdata.txt'), columns=['long', 'lat', 'energyb'])
dfa = pd.DataFrame(np.loadtxt('uvadata.txt'), columns=['long', 'lat', 'energya'])
temp = pd.DataFrame(np.loadtxt('daytemp.txt')) #temperature at Gale for solar longitude 0-360
daytemp = np.array(temp[1])
energy_c = np.array(dfc['energyc'])
energy_b = np.array(dfb['energyb'])
energy_a = np.array(dfa['energya'])
lat = np.array(dfc['lat'])
long = np.array(dfc['long'])
e_new_c = [None]*36
e_totalc = 0
e_new_b = [None]*36
e_totalb = 0
e_new_a = [None]*36
e_totala = 0
j = 0
Gale_lat = -4.5
#calculate UV energies throughout the year for Gale using interpolation
for i in range(len(energy_c)):
    if lat[i] == -10:
        x = [lat[i],lat[i+1]]
        y = [energy_c[i],energy_c[i+1]]
        e_new_c[j] = np.interp(Gale_lat,x,y)
        e_totalc = e_totalc + e_new_c[j]
        j = j+1
j = 0
for i in range(len(energy_b)):
    if lat[i] == -10:
        x = [lat[i],lat[i+1]]
        y = [energy_b[i],energy_b[i+1]]
        e_new_b[j] = np.interp(Gale_lat,x,y)
        e_totalb = e_totalb + e_new_b[j]
        j = j+1
j = 0
for i in range(len(energy_a)):
    if lat[i] == -10:
        x = [lat[i],lat[i+1]]
        y = [energy_a[i],energy_a[i+1]]
        e_new_a[j] = np.interp(Gale_lat,x,y)
        e_totala = e_totala + e_new_a[j]
        j = j+1
Ls = np.array(list(range(0,360,10))) #Ls = solar longitude


# In[19]:


years = np.array(np.linspace(1,50,num=50)) #time of 50 martian years
ddep = 3.5*10**-6 #3.5 micron yearly dust deposition at Gale
uva_abs = 29.35 #UV a absorption coefficient in cm-1 from Godin et al., 2023
uvb_abs = 29.54 #UV b absorption coefficient in cm-1
uvc_abs = 30.49 #UV c absorption coefficient in cm-1
energy_a = e_new_a
energy_b = e_new_b
energy_c = e_new_c
energy_time_c = np.array(np.zeros((len(years),len(Ls))))
energy_time_b = np.array(np.zeros((len(years),len(Ls))))
energy_time_a = np.array(np.zeros((len(years),len(Ls))))
IDPlife = np.array(np.zeros((len(years),len(IDPdiam)))) #IDP lifetime
etotal_Ls = np.array(np.zeros((len(years),len(Ls))))
ec = np.array([0]*len(years))
eb = np.array([0]*len(years))
ea = np.array([0]*len(years))
ec[0] = sum(energy_c)
eb[0] = sum(energy_b)
ea[0] = sum(energy_a)
energy_time_c[0,:] = energy_c
energy_time_b[0,:] = energy_b
energy_time_a[0,:] = energy_a
#calculate UV fluxes with time and gradual dust deposition (ddep)
for i in range(len(years)):
    for j in range(len(Ls)):
        if i>=0:
            energy_time_c[i,j] = energy_c[j]/(math.e**(i*uvc_abs*1000*ddep))
            ec[i] = ec[i] + energy_time_c[i,j]
            energy_time_b[i,j] = energy_c[j]/(math.e**(i*uvb_abs*1000*ddep))
            eb[i] = eb[i] + energy_time_b[i,j]
            energy_time_a[i,j] = energy_c[j]/(math.e**(i*uva_abs*1000*ddep))
            ea[i] = ea[i] + energy_time_a[i,j]


# In[20]:


'''Calculations for QE, IDP life, equations from Moores and Schuerger 2012 and Moores et al., 2017
QE(x,T)=4.87×10^(-14) T-2.73×10^(-12)
L_n=(2χD_n ρf)/(3M_w F_(UV, AVG) QE)'''

FUV_avg = 59354500/36 #averaged UV flux
seconds_Myr = 88775
X = 0.1 #10% Carbon weight percentage for IDPs
QE = np.array([None]*len(Ls))
FQEdaily = np.array(np.zeros((len(years),len(Ls))))
totalfqe = np.array([None]*len(years))
totaluv = np.array(np.zeros((len(years))))
IDPlife = np.array(np.zeros((len(years),len(IDPdiam))))
for i in range(len(Ls)):
    QE[i] = (daytemp[i]*((2.76*10**-15))-(1.54*10**-13)) 
for i in range(len(years)):
    for j in range(len(Ls)):
        totaluv[i] = (energy_time_c[i,j]+energy_time_b[i,j]+energy_time_a[i,j])
        FQEdaily[i,j] = QE[j]*(energy_time_c[i,j]+energy_time_b[i,j]+energy_time_a[i,j])
print(totaluv)
for i in range(len(years)):
    totalfqe[i] = 0
    for j in range(len(Ls)):
        totalfqe[i] = totalfqe[i]+FQEdaily[i,j]

totalfqe = totalfqe*FUV_avg/(seconds_Myr) 
for i in range(len(years)):
    for j in range(len(IDPdiam)):
        IDPlife[i,j] = (2*X*IDPunalt[j]*1000*IDPunmelt[j]*IDPdiam[j])/(3*totalfqe[i]*0.012) 


# In[21]:


carbon = np.array(np.zeros((len(years),len(IDPdiam))))
carbon[0,:] = (np.array(FlynnIN))
carbon[0,:] = carbon[0,:]
sumcarbon = np.array([None]*len(years))
layermass = np.zeros(len(years))

ppb = np.zeros(len(years))
frac = np.zeros(len(years))
for i in range(len(years)):
    layermass[i] = ddep*870 #density of JSC-1 Mars simulant 0.87 kg/cm3 (Allen et al., 1998)
    
for i in range(len(years)):
    for j in range(len(IDPdiam)):
        if i>=1:
            carbon[i,j] = (carbon[i-1,j]-carbon[i-1,j]/(IDPlife[i-1,j]))
for i in range(len(years)):
    sumcarbon[i] = sum(carbon[i,:])
    ppb[i] = (sumcarbon[i]/layermass[i])*10**9
for i in range(len(years)):
    frac[i] = (ppb[0]-ppb[i])/ppb[i] #fraction of organic carbon destroyed
plt.rc("font", family="serif", size=18.)
plt.figure(figsize =(9.5, 6.0))
plt.plot(years,frac*100,linestyle='-',linewidth=1.5,marker="d",color='darkred', markerfacecolor='#E69F00', markersize=5)
plt.xlabel('Martian Years',size=20)
plt.yscale('log')
plt.ylabel('% of organic carbon destroyed',size=20)
plt.title('Destruction of organics at Gale Crater',size=20)
plt.ylim(0.0,0.4)
plt.yticks([0.1,0.3], ['0.1','0.3'])
plt.tight_layout
plt.savefig('GaleOrganics_Das.png')


# In[ ]:





# In[ ]:




