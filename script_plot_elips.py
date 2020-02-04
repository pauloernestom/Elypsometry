#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 12:53:40 2018

@author: paulo
"""
import time
start_time = time.time()
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import os





path = #data_directory

plotDir = path + '/plots/'
pathTab = path + '/data/'

if not os.path.exists(plotDir):
    os.makedirs(plotDir)
if not os.path.exists(pathTab):
    os.makedirs(pathTab)

files = []
for i in os.listdir(path):
    if i.endswith(".dat"):
         files.append(path + str(i))
files.sort()

polymers_dict = {}
T_dict = {'T_80K':[], 'T_120K':['GO_RGO_150K'], 'T_150K':[], 'T_180K':[], 'T_210K':[], 'T_240K':[], 'T_270K':[], 'T_300K':[]}
#T_dict = {}
cols = ['Wavelength / nm']
for i in range(0, len(files)):
    cols.append(files[i].split('/')[-1][:-9])
    if files[i].split('/')[-1][:-9].rsplit('_',1)[0] in files[i]: 
        polymers_dict[files[i].split('/')[-1][:-9].rsplit('_',1)[0]] = list()
    
    
for i in range(0, len(files)):
    for b in polymers_dict:
        if files[i].split('/')[-1][:-9].rsplit('_',1)[0] == b:
            polymers_dict[files[i].split('/')[-1][:-9].rsplit('_',1)[0]].append(files[i].split('/')[-1][:-9])
        
    if (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) <= 100:
        T_dict['T_80K'].append(files[i].split('/')[-1][:-9])
        
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 120:
        T_dict['T_120K'].append(files[i].split('/')[-1][:-9])
        
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 150:
        T_dict['T_150K'].append(files[i].split('/')[-1][:-9])
    
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 180:
        T_dict['T_180K'].append(files[i].split('/')[-1][:-9])
        
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 210:
        T_dict['T_210K'].append(files[i].split('/')[-1][:-9])
        
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 240:
        T_dict['T_240K'].append(files[i].split('/')[-1][:-9])
        
    elif (int(files[i].split('/')[-1][:-9].rsplit('_',1)[1][:-1])) == 270:
        T_dict['T_270K'].append(files[i].split('/')[-1][:-9])
        
    else:
        T_dict['T_300K'].append(files[i].split('/')[-1][:-9])


cols_1=['Wavelength / nm']
for i in cols:
    if (i.find("PEDOT_PCBM")==0):
        cols_1.append(i)



cols_2=['Wavelength / nm_1']
for i in cols:
    if (i.find("PEDOT_PCBM")==-1) and (i.find("Wavelength")==-1):
        cols_2.append(i)

s1s0 = []

s1s0_1=[]
      

for line in files:
    if "PEDOT_PCBM" in line:
        s1s0_1.append(np.loadtxt(line, skiprows=1, usecols=(1,), unpack=True))
        break


s1s0.append(np.loadtxt(files[1], skiprows=1, usecols=(1,), unpack=True))


for line in range(0,len(files)):
    if "PEDOT_PCBM" in files[line]:
        s1s0_1.append(np.loadtxt(files[line], skiprows=1, usecols=(10,), unpack=True)[0:3631])
    else:
        s1s0.append(np.loadtxt(files[line], skiprows=1, usecols=(10,), unpack=True)[0:2048])

tab_s1s0_1 = pd.DataFrame(np.transpose(s1s0_1), columns=cols_1)
tab_s1s0 = pd.DataFrame(np.transpose(s1s0), columns=cols_2)


for i in cols_2:
    tab_s1s0_1[i]=tab_s1s0[i]


print((cols))

temperatures = ['T_80K', 'T_120K', 'T_150K', 'T_150K', 'T_180K', 'T_210K', 'T_270K', 'T_300K']

for i in polymers_dict:
    fig = plt.figure(figsize=(10,9))
    for l in polymers_dict[i]:
        if "PEDOT_PCBM" in i:
            plt.plot(tab_s1s0_1['Wavelength / nm'], tab_s1s0_1[l],  label=l)
        else:
            plt.plot(tab_s1s0_1['Wavelength / nm_1'], tab_s1s0_1[l],  label=l)
    plt.xlim(600,800)
    if i == 'GO_PCBM':
        plt.ylim(-.75,.75)
    elif i == 'GO_PCBM_RGO':
        plt.ylim(-1,1)
    elif i == 'GO_RGO':
        plt.ylim(-.5,.5)
    elif i =='PEDOT_GO_PCBM':
        plt.ylim(-1.5,1.5)
    elif i == 'PEDOT_GO_PCBM_RGO':
        plt.ylim(-1.2,1.2)
    elif i == 'PEDOT_GO_RGO':
        plt.ylim(-0.5,0.5)
    elif i == 'PEDOT_RGO':
        plt.ylim(-0.5,0.5)
    else:
        plt.ylim(-2,2)
    plt.xlabel('Wavelength / nm', size=16)
    plt.ylabel('Intensity / a.u.', size=16)
    plt.tick_params(direction='in', which='both', length=5)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    
    plt.legend(fontsize = 16)
#    plt.savefig(plotDir + 'fig-' + i + '-s1s0.eps', format='eps')
    plt.show()


for i in temperatures:
    fig = plt.figure(figsize=(10,9))
    for l in T_dict[i]:
        if "PEDOT_PCBM" in l:
            plt.plot(tab_s1s0_1['Wavelength / nm'], tab_s1s0_1[l],  label=l)
        else:
            plt.plot(tab_s1s0_1['Wavelength / nm_1'], tab_s1s0_1[l],  label=l)
    plt.xlim(600,800)
    if i == 'T_80K':
        plt.ylim(-1.5,1.5)
    elif i == 'T_120K':
        plt.ylim(-2,1)
    elif i == 'T_150K':
        plt.ylim(-1.5,.5)
    elif i == 'T_210K':
        plt.ylim(-1.2,1.2)
    elif i == 'T_240K':
        plt.ylim(-1,1)
    else:
        plt.ylim(-2,2)
    plt.xlabel('Wavelength / nm', size=16)
    plt.ylabel('Intensity / a.u.', size=16)
    plt.tick_params(direction='in', which='both', length=5)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    
    plt.legend(fontsize = 16)
#    plt.savefig(plotDir + 'fig-' + i + '-s1s0.eps', format='eps')
    plt.show()


#save data 


"""   

for i in polymers:
    if "PEDOT_PCBM" in i:
        f2 = open(pathTab + i + '_s1s0' + '.dat', 'w')
        tab_s1s0_1[['Wavelength / nm']+ globals()[i]].copy().to_csv(f2, sep='\t')
        f2.close()
    else:
        f2 = open(pathTab + i + '_s1s0' + '.dat', 'w')
        tab_s1s0_1[['Wavelength / nm_1']+ globals()[i]].copy().to_csv(f2, sep='\t')
        f2.close()

for i in temperatures:
    f2 = open(pathTab + i + '_s1s0' + '.dat', 'w')
    tab_s1s0_1[['Wavelength / nm_1'] + ['Wavelength / nm'] + globals()[i]].copy().to_csv(f2, sep='\t')
    f2.close()
"""
print('\n  --- %s seconds --- \n'% (time.time() - start_time))



########