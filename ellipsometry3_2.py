# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 17:19:46 2018

@author: Paulo
"""
import time
start_time = time.time()
import os
import sys
import pandas as pd
import Ellipsometry_calc as elc






path1 = os.getcwd()
path = str(sys.argv[1:])[2:-2] #read the path from terminal


######################################################################
#crate a list of the directories
######################################################################
dir_list = os.listdir(path)

listdir=[]
for d in dir_list:
    listdir.append(path + d)

tb_lisdir = pd.DataFrame(data=listdir)

f1 = open('tb_lisdir.txt', 'w')
tb_lisdir[0].to_csv(f1, sep='\t')
f1.close()

#
print('path:', path)

_filext = '.txt'



for c in range(0,len(listdir)):
    print(listdir[c])
    figname = str(listdir[c].split('/')[-1] + '_2')

    dataname = str(listdir[c].split('/')[-1] + '_data')

   
    rd = elc.Read_Data(filename=_filext ,path=str(listdir[c]+'/'))


    files=rd.loadData()
    data=rd.readData_elipsometry()
    wavelength=rd.readData_elipsometry(wavelength=True)
    
    _num_of_files = len(files)
    
    
    
    pathTab=rd.dirdata()
    plotDir=rd.dirplot()

    degrees = rd.degree_read()

    calc = elc.Manage_Data(data=data, degrees=degrees)
    
    A=calc.Fourier_terms('A')
    
    B=calc.Fourier_terms('B')
    
    C=calc.Fourier_terms('C')
    
    D=calc.Fourier_terms('D')
    
    
    
    S0 = calc.Stokes_parameters('S0')

    S1 = calc.Stokes_parameters('S1')
    
    S2 = calc.Stokes_parameters('S2')
    
    S3 = calc.Stokes_parameters('S3')

    s1s0 = calc.Div_ss('s1s0')
    
    s2s0 = calc.Div_ss('s2s0')
    
    s3s0 = calc.Div_ss('s3s0')

    P = calc.Other_facs('P')

    r = calc.Other_facs('r')
    
    g = calc.Other_facs('g')
    
    Phi = calc.Other_facs('Phi')
    
    Chi = calc.Other_facs('Chi')
   

    
    sf = elc.Save_Data(wavelength, A, B, C, D,S0, S1,S2, S3, s1s0,s2s0,s3s0, P, r, g, Phi, Chi)
    

    sf.save_tab(pathTab, dataname)

    ptdta = elc.Plot_Data(wavelength, A, B, C, D,S0, S1,S2, S3, s1s0,s2s0,s3s0, P, r, g, Phi, Chi)
    
    ptdta.plotes('stokes', plotDir, figname, save_fig=True)
    ptdta.plotes('factors', plotDir, figname, save_fig=True)

#for c in range(0,len(listdir)):
#    print(listdir[c])
#    figname = str(listdir[c].split('/')[-1]+ '_temp')



   

print('\n  --- %s seconds --- \n'% (time.time() - start_time))
