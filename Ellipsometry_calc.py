# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:28:17 2018

@author: Paulo
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pylab as plt



class Read_Data():
    """docstring for ReadData """
    def __init__(self, filename='', path='', info=None):
        self.filenam = filename
        self.path = path
        self.info = info
        self.data=[]
        self.files = []
        self.files_info = []
        self.cols = []
		
    def dirplot(self):
        self.plotDir = self.path + 'plots/'
        if not os.path.exists(self.plotDir):
            os.makedirs(self.plotDir)
        return self.plotDir
    def dirdata(self):
        self.pathTab = self.path + 'data/'
        if not os.path.exists(self.pathTab):
            os.makedirs(self.pathTab)
        return self.pathTab
        
            
    def loadData(self):

        for i in os.listdir(self.path):
            if i.endswith(self.filenam):
                self.files.append(self.path + str(i))
        self.files.sort(key=lambda x: int(x.split('/')[-1][0:-4]))
       
        return self.files
    
    def readData_elipsometry(self, wavelength=False):
        if wavelength:
            for i in self.files:
                self.data = (np.loadtxt(i, comments='>', skiprows=17, usecols=(0,), unpack=True))
                break
            
        else:
            for i in self.files:
                self.data.append(np.loadtxt(i, comments='>', skiprows=17, usecols=(1,), unpack=True))
			
        return self.data	
    
    def degree_read(self):
        self.degree=[]
        for i in self.files:
            self.degree.append(float(i.split('/')[-1][0:-4]))
        return self.degree

class Manage_Data():
    """docstring for Manage_Data """
    def __init__(self, data='', degrees=''):
        self.data=data
        self.degrees = degrees
        
    def Fourier_terms(self, term): #Valid term: A, B, C or D
        if term == 'A':
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _Term=sum(map(A_fuc,self.data))
        elif term == 'B':
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _Term=sum(map(B_fuc,self.data,self.degrees))
        elif term == 'C':
            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _Term=sum(map(C_fuc,self.data,self.degrees))
        elif term == 'D':
            def D_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(4*deg_i *(np.pi/180)))
            _Term=sum(map(D_fuc,self.data,self.degrees))
        else:
            _Term=print(str(term)+' is not a valid term! Try A, B, C or D')
        return _Term
    
    def Stokes_parameters(self, parm): #Valid parm: S0, S1, S2 or S3
        
        if parm == 'S0':
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _Parm = _A - _C
        
        elif parm == 'S1':
          
            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            

            _Parm = -2*_C
        
        elif parm == 'S2':
            
            def D_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(4*deg_i *(np.pi/180)))
            _D=sum(map(D_fuc,self.data,self.degrees))
            
            _Parm = -2*_D
        
        elif parm == 'S3':
         
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _B=sum(map(B_fuc,self.data,self.degrees))

            _Parm = -_B
        else:
            _Parm=print(str(parm)+' is not a valid term! Try S0, S1, S2 or S3')
        
        return _Parm
    
    def Div_ss(self, div): # Valid div: s1s0, s2s0 or s3s0
        
        def div_ss(S, s):
            def div_check(x, y):
              try:
                x / y
              except ZeroDivisionError:
                return 0
              else:
                return x/y
            x_y=[]    
    
            for (i, f) in zip(S, s):
                x_y.append(div_check(float(i),float(f)))
            return x_y
        
        if div == 's1s0':            
            
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C
            
            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            

            _S1 = -2*_C
            
            _Div = div_ss(_S1, _S0)
        
        elif div == 's2s0':
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C
            
            
            def D_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(4*deg_i *(np.pi/180)))
            _D=sum(map(D_fuc,self.data,self.degrees))
            
            _S2 = -2*_D

            _Div = div_ss(_S2, _S0)
            
        elif div == 's3s0':
            
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _B=sum(map(B_fuc,self.data,self.degrees))

            _S3 = -_B
            
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C

            _Div = div_ss(_S3, _S0)

        else:
            _Div=print(str(div)+' is not a valid term! Try s1s0, s2s0 or s3s0.')
        
        return _Div
        
#Other factors
        #Degree of polarization - P
        #Anisotropy - r
        #Asymmetry - g 
        #Angle Phi and ellipticity Chi
        
    def Other_facs(self, fac): #Valid fac: P, r, g, Phi or Chi
        def power(my_list):
            return [ x**2 for x in my_list ]
        if fac == 'P':
            def div_ss(S, s):
                def div_check(x, y):
                  try:
                    x / y
                  except ZeroDivisionError:
                    return 0
                  else:
                    return x/y
                x_y=[]    
                for i, f in zip(S, s):
                    x_y.append(div_check(float(i),float(f)))
                return x_y
            
                
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C
            
            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            

            _S1 = -2*_C
            
            s1s0 = div_ss(_S1, _S0)
        

            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C
            
            
            def D_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(4*deg_i *(np.pi/180)))
            _D=sum(map(D_fuc,self.data,self.degrees))
            
            _S2 = -2*_D

            s2s0 = div_ss(_S2, _S0)
            
            
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _B=sum(map(B_fuc,self.data,self.degrees))

            _S3 = -_B
            
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C

            s3s0 = div_ss(_S3, _S0)
#
#            def power(my_list):
#                return [ x**2 for x in my_list ]
            Fac = np.sqrt(power(s1s0)) + (power(s2s0)) + (power(s3s0))
        
        elif fac == 'r':
            def div_ss(S, s):
                def div_check(x, y):
                  try:
                    x / y
                  except ZeroDivisionError:
                    return 0
                  else:
                    return x/y
                x_y=[]    
                for i, f in zip(S, s):
                    x_y.append(div_check(float(i),float(f)))
                return x_y
            
                
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C
            
            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            

            _S1 = -2*_C
            
            s1s0 = div_ss(_S1, _S0)
            
            s1s02=[]
            s1s23=[]
            for i in s1s0:
                s1s02.append(-2.0*i)
            for i in s1s0:
                s1s23.append(3.0+i)
            Fac=[]
            for i,f in zip(s1s02,s1s23):
                Fac.append(i/f)
            #Fac = (s1s02)/(s1s23)
            
        elif fac == 'g':
            def div_ss(S, s):
                def div_check(x, y):
                  try:
                    x / y
                  except ZeroDivisionError:
                    return 0
                  else:
                    return x/y
                x_y=[]    
                for i, f in zip(S, s):
                    x_y.append(div_check(float(i),float(f)))
                return x_y
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _B=sum(map(B_fuc,self.data,self.degrees))

            _S3 = -_B
            
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C

            s3s0 = div_ss(_S3, _S0)
            
            Fac=[]
            for i in s3s0:
                Fac.append(2.0*i)

 

#to aqui!       
        elif fac == 'Phi':
            def div_ss(S, s):
                def div_check(x, y):
                  try:
                    x / y
                  except ZeroDivisionError:
                    return 0
                  else:
                    return x/y
                x_y=[]    
                for i, f in zip(S, s):
                    x_y.append(div_check(float(i),float(f)))
                return x_y
            def D_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(4*deg_i *(np.pi/180)))
            _D=sum(map(D_fuc,self.data,self.degrees))
            
            _S2 = -2*_D

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            

            _S1 = -2*_C
            
            s2s1 = div_ss(_S2, _S1)
            
            Fac = np.arctan(s2s1)

        elif fac == 'Chi':
            def div_ss(S, s):
                def div_check(x, y):
                  try:
                    x / y
                  except ZeroDivisionError:
                    return 0
                  else:
                    return x/y
                x_y=[]    
                for i, f in zip(S, s):
                    x_y.append(div_check(float(i),float(f)))
                return x_y
            def B_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.sin(2*deg_i *(np.pi/180)))
            _B=sum(map(B_fuc,self.data,self.degrees))

            _S3 = -_B
            
            def A_fuc(data_i):
                return (2/9)*(data_i)
            _A=sum(map(A_fuc,self.data))

            def C_fuc(data_i, deg_i):
                return (4/9)*(data_i*np.cos(4*deg_i *(np.pi/180)))
            _C=sum(map(C_fuc,self.data,self.degrees))            
            
            _S0 = _A - _C

            s3s0 = div_ss(_S3, _S0)
            
            Fac = np.arctan(s3s0)
        
        else:
            Fac=print(str(fac)+' is not a valid term! Try P, r, g, Phi or Chi.')
            
        return Fac
        

    
class Save_Data():
    """docstring for Save_Data """
    def __init__(self, wavelength, A, B, C, D,S0, S1,S2, S3, s1s0,s2s0,s3s0, P, r, g, Phi, Chi):
        self.wavelength = wavelength
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.S0 = S0
        self.S1 = S1
        self.S2 = S2
        self.S3 = S3
        self.s1s0 = s1s0
        self.s2s0 = s2s0
        self.s3s0 = s3s0
        self.P = P
        self.r = r
        self.g= g
        self.Phi = Phi
        self.Chi = Chi
        self.data = []
        self.tab_data = []

    def save_tab(self, directory, dataname):
        self.cols=['wavelength', 'A', 'B', 'C', 'D', 'S0', 'S1', 'S2', 'S3', 's1s0', 's2s0', 's3s0', 'P', 'r', 'g', 'Phi', 'Chi']
        self.data.append(self.wavelength)
        self.data.append(self.A)
        self.data.append(self.B)
        self.data.append(self.C)
        self.data.append(self.D)
        self.data.append(self.S0)
        self.data.append(self.S1)
        self.data.append(self.S2)
        self.data.append(self.S3)
        self.data.append(self.s1s0)
        self.data.append(self.s2s0)
        self.data.append(self.s3s0)
        self.data.append(self.P)
        self.data.append(self.r)
        self.data.append(self.g)
        self.data.append(self.Phi)
        self.data.append(self.Chi)
        self.tab_data=pd.DataFrame(data=np.transpose(self.data), columns=self.cols)

   
        self.directory = directory
        self.dataname = dataname 
        self.f2 = open(self.directory + self.dataname + '.dat', 'w')
        self.tab_data.to_csv(self.f2, sep='\t')
        self.f2.close()
    


        
class Plot_Data():
    """docstring for Plot_Data """
    def __init__(self, wavelength, A, B, C, D,S0, S1,S2, S3, s1s0,s2s0,s3s0, P, r, g, Phi, Chi):
        self.wavelength = wavelength
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.S0 = S0
        self.S1 = S1
        self.S2 = S2
        self.S3 = S3
        self.s1s0 = s1s0
        self.s2s0 = s2s0
        self.s3s0 = s3s0
        self.P = P
        self.r = r
        self.g= g
        self.Phi = Phi
        self.Chi = Chi
        self.data = []
        self.tab_data = []

    def plotes(self, plot, directory, figname, save_fig=False, show=False):#valid:stokes, factors
        self.show = show
        self.plot = plot
        self.directory = directory
        self.figname = figname
        self.save_fig = save_fig
        def plot_config(xlim, ylim, xlabel, ylabel):
            plt.legend(fontsize=16)
            plt.xlim(xlim)#550,700
            plt.xticks(fontsize=16)
            plt.tick_params(direction='in', which='major', length=7)
            plt.ylim(ylim)#-1.5,1
            plt.yticks(fontsize=16)
            plt.ylabel(xlabel, size=16)#'Wavelength / nm'
            plt.xlabel(ylabel, size=16)#'Intensity / a.u.'
        
        if self.plot == 'stokes':
            self.fig=plt.figure(figsize=(10, 8))
            plt.plot(self.wavelength,self.s1s0, label='S$_{1}$/S$_{0}$')
            plt.plot(self.wavelength,self.s2s0, label='S$_{2}$/S$_{0}$')
            plt.plot(self.wavelength,self.s3s0, label='S$_{3}$/S$_{0}$')
            plot_config([600,800],[-1.5,1],'Intensity / a.u.','Wavelength / nm')
            if self.show:
                plt.show()
            if self.save_fig:
                plt.savefig(self.directory + self.figname + '.eps', format='eps')
        
        elif self.plot == 'factors':
            self.fig=plt.figure(figsize=(10, 8))
            plt.plot(self.wavelength,self.P, label='Degree of polarization (P)')
            plt.plot(self.wavelength,self.r, label='Anisotropy (r)')
            plt.plot(self.wavelength,self.g, label='Asymmetry (g)')
            plot_config([600,800],[-1,1.5],'Intensity / a.u.','Wavelength / nm')
            if self.show:
                plt.show()
            if self.save_fig:
                plt.savefig(self.directory + self.figname + '_factors' + '.eps', format='eps')  
                                  
        else: 
            print(str(plot)+' is not a valid term! Try stokes or factors')
            
    def plot_temp(self, s3s0, directory, figname, label, save_fig=False, show=False):
        def plot_config(xlim, ylim, xlabel, ylabel):
            plt.legend(fontsize=16)
            plt.xlim(xlim)#550,700
            plt.xticks(fontsize=16)
            plt.tick_params(direction='in', which='major', length=7)
            plt.ylim(ylim)#-1.5,1
            plt.yticks(fontsize=16)
            plt.ylabel(xlabel, size=16)#'Wavelength / nm'
            plt.xlabel(ylabel, size=16)#'Intensity / a.u.'
        self.fig=plt.figure(figsize=(10, 8))
        plt.plot(self.wavelength,self.s3s0, label=label)
        plot_config([600,800],[-1.5,1],'Intensity / a.u.','Wavelength / nm')
        if self.show:
            plt.show()
        if self.save_fig:
            plt.savefig(self.directory + self.figname + '.eps', format='eps')