# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:57:58 2016

@author: haghakha
"""

import numpy as np
import os
import shutil as st
import pyDOE as pylhs
import subprocess as sprcs
import time


num_sample=64
randoms=2

lhd = pylhs.lhs(randoms, num_sample)
lhs=np.zeros((num_sample,randoms))
lhs[:,0]=lhd[:,0]*.04+.01
lhs[:,1]=lhd[:,1]*400+800

np.savetxt("lhs_samples.csv", lhs, delimiter=",")

fp = open(os.path.join(os.getcwd(),'puffin/puffin.inp'),'r')
simdata = fp.read()
fp.close


water_w=0.33
temp=1111

dircs_name='LHS'
#delete the folders if they are already exist
for i in range(0,num_sample): 
    dirname=os.path.join(os.getcwd(),dircs_name+str(i))
    if (os.path.isdir(dirname)):
        st.rmtree(os.path.join(os.getcwd(),dircs_name+str(i))) 

#creates the folder for simulation
for i in range(0,num_sample): 
    st.copytree('puffin',os.path.join(os.getcwd(),dircs_name+str(i))) 

#replace the sample value with initial value in the list of input    
for i in range(0,num_sample):        
    print lhs[i,:]
    rep_water_w="{:.4f}".format(lhs[i,0])
    rp_simdata = simdata.replace(str(water_w),rep_water_w)
    rep_temp="{:.0f}".format(lhs[i,1])
    rp_simdata = rp_simdata.replace(str(temp), rep_temp)      
    fp = open(os.path.join(os.getcwd(),dircs_name+str(i)+'/puffin.inp'),'w')
    fp.write(rp_simdata)
    fp.flush()
    fp.close
    print i

for i in range(0,num_sample):
    path=dircs_name+str(i)
    os.chdir(path)    
    os.system('./puffin>output')
    os.system('gnuplot field.gnu')
    image = sprcs.Popen(["display", "./conc.jpg"])
    time.sleep(1)
    image.kill()    
    os.chdir('..')
    print i
    
particle_flux=np.zeros(num_sample)
eruption_height=np.zeros(num_sample)
    
for i in range(0,num_sample): 
    dirname=os.path.join(os.getcwd(),dircs_name+str(i))
    fp = open(os.path.join(dirname,'output'),'r')
    for line in fp:
        if line.find("PARTICLE FLUX AT HB")>-1:
            info = line.split()
            particle_flux[i] = float(info[4])
        if line.find("ERUPTION PLUME HEIGHT")>-1:
            info = line.split()
            eruption_height[i]=float(info[3]) 
            
h_mean=sum(eruption_height)/num_sample  
h_std=np.sqrt(sum(np.square(eruption_height))/num_sample-h_mean*h_mean)  

np.savetxt("eruption_height.csv", eruption_height, delimiter=",") 
np.savetxt("particle_flux.csv", particle_flux, delimiter=",")


print eruption_height
print h_mean
print h_std
                         
           
            
    
    
    