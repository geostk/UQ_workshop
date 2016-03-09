# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:22:43 2016

@author: zhixuanc
"""

import numpy as np
import scipy.special as sps
import os
import shutil as st
import subprocess as sprcs
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
"""
==============================================================================================SUB FUNCTIONS: FUNCTIONS THAT WILL BE CALLED IN MAIN FUNCTION============================================================================================================

===============================================================================================================DO NOT MODIFY THIS CODE============================================================================================================
"""
"""
 Recurrence coefficients for monic Jacobi polynomials.
 This function will be called to generate 
"""
def jacobi(N,a,b):

    nu=(b-a)/(a+b+2.0)
    mu=2**(a+b+1)*sps.gamma(a+1)*sps.gamma(b+1)/sps.gamma(a+b+2)
    if N==1:
        ab=np.append(nu, mu)
    else:
        N=N-1
        n=np.array(range(1,N+1))
        nab=2*n+float(a)+b
        A=np.append(nu, np.divide((b**2-a**2)*np.ones(N), np.multiply(nab, (nab+2))))
        
        B1=4*(float(a)+1)*(b+1)/((a+b+2)**2*(a+b+3))
        C=np.append(mu,B1)
        if N>1:
            n=np.array(range(1,N))
            nab=nab[n]
            n+=1
            mul1=np.multiply((n+float(a)),(n+b))
            mul2=np.multiply(n, (n+float(a)+b))
            mul3=np.multiply(mul1, mul2)
            mul4=np.multiply(np.multiply((nab**2),(nab+1)),(nab-1))
            B=4*np.divide(mul3,mul4)
            D=np.append(C,B)
            E=np.array([A,D])        
        else:
            E=np.array([A,C])
        ab=E.T
    return ab 
    
"""
function that should return Ph(si, j) for Beta distribution
 si is a series of points in random space 
 j indicate the order of polynomials 
"""

def jacobi_poly (si,c,alfa, beta):
    phi=[0.]*len(si)
    for i in range(0, len(c)):
        if c[i] != 0:
            for j in range(0, len(si)): 
                phi[j] +=c[i]*sps.eval_jacobi(i, alfa, beta, si[j])
            
    return np.array(phi) 
    
"""
 function that should return Ph(si, j) for XX distribution
 si is a series of points in random space, should be a list, array or tuple
 j indicate the order of polynomials 
 return values of j th order polynomials at points si 
 
 Note: assuming alf and beta to be zero for Jacobi polynomial
"""   
def XX_poly(si,j,dis, alfa, beta):
    c=[0]*len(si)
    
    for i in range(0,len(si)):
        if i==j:
            c[i]=1
        else:
            c[i]=0    
    
    if dis=="Gaussian":
        phi=np.polynomial.hermite_e.hermeval(si, c)
    elif dis=="Gamma":
        phi=np.polynomial.laguerre.lagval(si, c)
    elif dis=="Beta":
        phi=jacobi_poly (si,c,alfa, beta)    
    elif dis=="Uniform":
        phi=np.polynomial.legendre.legval(si, c)
    else:
         print ("input parameter dis(distribution of underlying random variables) is not recognizable, please double check!")   
    
    
    return phi

"""
 function that determine the coefficients based on PCQ
 h should be 2 dimension array m by n: m is number of samples in 1 direction, n is number of samples in another direction  
 w should be 2 dimension array m by n: m is number of samples in 1 direction, n is number of samples in another direction
 si is a list of quadrature points in each dimension
 dis is the name of ditribution in that direction
 """
def calcu_coefficient (dis1,dis2,h,w,si1,si2, alf, beta, v_name):
    basis=np.zeros((len(si1),len(si2)))
    coef= np.array([[0]*len(si2)]*len(si1))
    for m in range(0,len(si1)):
        Phi=np.mat(XX_poly(si1,m,dis1,alf,beta))
        for n in range(0,len(si2)):
            Phj=np.mat(XX_poly(si2,n,dis2,alf,beta))
            basis=np.array(Phi.T*Phj)
            sum_num =sum(sum(np.multiply(np.multiply(h,basis),w)))
            sum_den =sum(sum(np.multiply(np.multiply(basis,basis),w)))     
            coef[m][n]=sum_num/sum_den
    
    print ("These coefficients for "+ v_name + " in front of basis functions are:")
    print coef
    
    return coef
    
"""
function for determining coefficient by solveing system of equations, Not PCQ at all
But more general
Will ber used when "Orthogonality" of basis can not been guaranteed
Questions: Is there any method that can give us "quadrature points, corresponding weight" and a set of orthogonal basis ----> will makes life easier 
"""
def calcu_coef_solve_Matrix(h, si1, si2, dis1, dis2, alf, beta, v_name):
    mat_size=len(si1)*len(si2)
    h_c=h.reshape((mat_size,1))
    A=np.mat(np.zeros((len(si1), len(si2))))
    B=np.mat(np.zeros((mat_size, 1)))
    M=np.mat(np.zeros((mat_size, mat_size)))
    
    count =0
    for i in range(0, len(si1)):
        Phi=np.mat(XX_poly(si1,i,dis1, alf, beta))
        for j in range(0, len(si2)):
            Phj=np.mat(XX_poly(si2,j,dis2, alf, beta))
            A=Phi.T*Phj
            B=A.reshape((len(si1)*len(si2),1))
            i
            j
            M[:, count]=np.array(B)
            count +=1
            
    coef=np.array(np.linalg.solve(np.array(M), h_c))
    
    print ("These coefficients for "+ v_name + " in front of basis functions are:")
    print coef
    
    return coef.reshape((len(si1), len(si2)))    
    
"""
function for evaluate height (or other output properties) 
input: dis, the name of the distribution
       si, the smaple points of underlying (standard) RVs
       coef, coefficients for each basis
"""
def test_coef(coef, si1, si2, dis1, dis2, alf, beta):
    h_out=np.zeros((len(si1),len(si2)))
    PHI=np.zeros((len(si1),len(coef),len(coef[0])))
    PHJ=np.zeros((len(si2),len(coef),len(coef[0])))
    
    for m in range(0,len(coef)):
        Phi=XX_poly(si1,m,dis1,alf, beta)
        for n in range(0,len(coef[0])):
            Phj=XX_poly(si2,n,dis2,alf, beta)
            PHI[:,m,n]=Phi
            PHJ[:,m,n]=Phj
            
    for i in range(0,len(si1)):
        for j in range(0,len(si2)):
             MI=(PHI[i,:,:])
             MJ=(PHJ[j,:,:])
             h_out[i,j]=sum(sum(np.multiply(np.multiply(MI,MJ),coef)))
    return h_out
    
"""
 Run simulation with sampled parameters
 input: sample_radius: sample points in first random direction
        sample_speed:  sample points in second random direction
        radius : value of first random parameter in original input file for puffin
        speed  : value of first random parameter in original input file for puffin
        dircs_name: The name of the method that you are using, for example 'PCQ'
"""
def run_puffin (samples_radius, samples_speed, radius, speed, dircs_name):
    num_sample_radius =len(samples_radius)
    num_sample_speed = len(samples_speed)
    total_samples = num_sample_radius*num_sample_speed
    
    fp = open(os.path.join(os.getcwd(),'puffin/puffin.inp'),'r')
    simdata = fp.read()
    fp.close

    #delete the folders if they are already exist
    for i in range(0,total_samples): 
        dirname=os.path.join(os.getcwd(),dircs_name+str(i))
        if (os.path.isdir(dirname)):
            st.rmtree(os.path.join(os.getcwd(),dircs_name+str(i))) 

    #creates the folder for simulation
    for i in range(0,total_samples): 
        st.copytree('puffin',os.path.join(os.getcwd(),dircs_name+str(i))) 

    #replace the sample value with initial value in the list of input    
    for i in range(0,num_sample_radius):      
        for j in range(0,num_sample_speed): 
            rep_radius="{:.0f}".format(samples_radius[i])
            rp_simdata = simdata.replace(str(radius),rep_radius)
            rep_speed="{:.0f}".format(samples_speed[j])
            rp_simdata = rp_simdata.replace(str(speed), rep_speed)      
            fp = open(os.path.join(os.getcwd(),dircs_name+str(i*num_sample_radius+j)+'/puffin.inp'),'w')
            fp.write(rp_simdata)
            fp.flush()
            fp.close   

    for i in range(0,total_samples):
        path=dircs_name+str(i)
        os.chdir(path)    
        os.system('./puffin>output')
        os.system('gnuplot field.gnu')
        image = sprcs.Popen(["display", "./conc.jpg"])
        time.sleep(1)
        image.kill()
        os.chdir('..')
        print i  

    print "Running of puffin is finished!"
    
"""
Parse output file and extract height and particle mass flux
input is: total_samples
will return (particle_flux , eruption_height)
"""    
def extract_height_particle_flux (total_samples, dircs_name):    
    
    particle_flux=np.zeros(total_samples)
    eruption_height=np.zeros(total_samples)
    
    for i in range(0,total_samples): 
        dirname=os.path.join(os.getcwd(),dircs_name+str(i))
        fp = open(os.path.join(dirname,'output'),'r')
        for line in fp:
            if line.find("PARTICLE FLUX AT HB")>-1:
                info = line.split()
                particle_flux[i] = float(info[4])
            if line.find("ERUPTION PLUME HEIGHT")>-1:
                info = line.split()
                eruption_height[i]=float(info[3])
    
    out_file_name1=dircs_name+"UQ_eruption_height.csv"
    out_file_name2=dircs_name+"UQ_particle_flux.csv"      
    np.savetxt(out_file_name1, eruption_height, delimiter=",") 
    np.savetxt(out_file_name2, particle_flux, delimiter=",")
                
    return (particle_flux, eruption_height)
                
    print "extract height and particle flux, done!"
 
"""
function that will generate sample points based on input distribution name
input: dis, distribution name
       num_sample: number of sample points
       alf, beta, parameters for beta distribution
"""   
def sample_gaussian_dis(dis, num_sample, alf, beta):
    if dis=="Gaussian":
        smpling=np.polynomial.hermite.hermgauss(num_sample)
    elif dis=="Gamma":
        smpling=np.polynomial.laguerre.laggauss(num_sample)
    elif dis=="Beta":
        smpling=jacobi(num_sample, alf, beta)    
    elif dis=="Uniform":
        smpling=np.polynomial.legendre.leggauss(num_sample)
    else:
         print ("input parameter dis(distribution of underlying random variables) is not recognizable, please double check!") 
    
    return smpling

"""
Function that will generate random samples for specific stabdard distribution
    parameters should keep consistent with the distribution that used to determine othogonal basis: the standard parameters
    Because these used to generate quadrature points should also be standard

"""
def generate_samples(dis, number_of_smaples_each_rv, a, b):
    if dis=="Gaussian":
        smpling=np.random.normal(0,1, number_of_smaples_each_rv)
    elif dis=="Gamma":
        smpling=np.random.standard_gamma(1, number_of_smaples_each_rv)
    elif dis=="Beta":
        smpling=np.random.beta(a,b, number_of_smaples_each_rv)   
    elif dis=="Uniform":
        smpling=np.random.uniform(-1, 1, number_of_smaples_each_rv)
    else:
         print ("input parameter dis(distribution of underlying random variables) is not recognizable, please double check!") 
    
    return smpling
"""
Wrapp up of PCQ
"""  
def PCQ_UQ (num_sample, mean_rv, range_rv, v_dis, name, rv_value, alf, beta):

    num_sample_radius=num_sample[0]
    num_sample_speed=num_sample[1]
    total_samples=num_sample_radius*num_sample_speed

    radius=rv_value[0]
    speed=rv_value[1]

    mean_radius=mean_rv[0]
    range_radius=range_rv[0]

    mean_speed=mean_rv[1]
    range_speed=range_rv[1]

    dis1=v_dis[0]
    dis2=v_dis[1]
    
    smpling1=sample_gaussian_dis(dis1, num_sample_radius, alf, beta)
    smpling2=sample_gaussian_dis(dis2, num_sample_speed, alf, beta)

    #This portition might not be correctly done
    samples_radius=range_radius*smpling1[0]+mean_radius
    samples_speed=range_speed*smpling2[0]+mean_speed

    weight1=smpling1[1]
    weight2=smpling2[1]

    """
    ==============================================================================================modification needed below this line============================================================================================================
    """
    dircs_name ='PCQ'
    run_puffin (samples_radius, samples_speed, radius, speed, dircs_name)

    particle_flux, eruption_height = extract_height_particle_flux (total_samples, dircs_name)
    """
    Compute mean and standard variance
    """ 
    """                                  
    h_mean=0
    h_sqr=0
    particle_flux_mean=0
    particle_flux_sqr=0
    for i in range(0,num_sample_radius):      
        for j in range(0,num_sample_speed): 
            h_mean=h_mean+weight1[i]*weight2[j]*eruption_height[i*num_sample_radius+j]
            h_sqr=h_sqr+weight1[i]*weight2[j]*eruption_height[i*num_sample_radius+j]*eruption_height[i*num_sample_radius+j]
            particle_flux_mean=particle_flux_mean+weight1[i]*weight2[j]*particle_flux[i*num_sample_radius+j]
            particle_flux_sqr=particle_flux_sqr+weight1[i]*weight2[j]*particle_flux[i*num_sample_radius+j]*particle_flux[i*num_sample_radius+j]

    h_std=np.sqrt(h_sqr-h_mean*h_mean) 
    particle_flux_std= np.sqrt(particle_flux_sqr-particle_flux_mean*particle_flux_mean)
        
    print 'vector of eruption height:', eruption_height
    print 'h_mean:', h_mean  
    print 'std dev for h:', h_std
    print 'vector of particle flux:',particle_flux
    print 'particle_flux_mean:  ',particle_flux_mean     
    print 'std dev for particle_flux:', particle_flux_std
    
    """

    """
    Prepare for coefficient computing
    """           
    h_new=eruption_height.reshape(num_sample_radius, num_sample_speed)
    p_new=particle_flux.reshape(num_sample_radius, num_sample_speed)
    si1=(smpling1[0]).tolist()
    si2=(smpling2[0]).tolist()

    w_new=np.zeros((num_sample_radius, num_sample_speed))
    for i in range(0,num_sample_radius):
        for j in range(0,num_sample_speed):
            w_new[i][j]=weight1[i]*weight2[j]

    """
    compute coefficient by PCQ
    The following way works perfectly for both RVs have the same distribution or the same number of points ---(just tested with 4 and 5 sample points, can not guarantee always works)
    Any way, the prerequirements that PCQ works for multi-dimension of RV is: the basis that generated based on tensor product of polynomial basis in each direction should be a set of orthogonal basis
    """
    coef=calcu_coefficient (dis1, dis2, h_new, w_new, si1, si2, alf, beta, "height")

    h_out=test_coef(coef, si1, si2, dis1, dis2, alf, beta)

    coef_p=calcu_coefficient (dis1, dis2, p_new, w_new, si1, si2, alf, beta, "particle flux")

    p_out=test_coef(coef_p, si1, si2, dis1, dis2, alf, beta)

    """
    A more general way to get coefficients 
    But not a PCQ at all
    Does not require orthogonality of basis
    """

    #coef=calcu_coef_solve_Matrix(h_new, si1, si2, dis1, dis2, alf, beta, "height")

    #h_out=test_coef(coef, si1, si2, dis1, dis2, alf, beta)

    #coef_p=calcu_coef_solve_Matrix(p_new, si1, si2, dis1, dis2, alf, beta, "particle flux")

    #p_out=test_coef(coef_p, si1, si2, dis1, dis2, alf, beta)

    print "height from running puffin"
    print h_new
    
    print "calculated height based on coefficient"
    print h_out

    print "particle mass flux from running puffin"
    print p_new
    
    print "calculated particle mass flux based on coefficient"
    print p_out
    """
    
    We plot h as a function of random parameters here:
    1) generate sample point ---> convert the smaple point in standard range into parameter range
    2) substitute these ample points into approximation expression based on calculated coefficient---> get corresponding height
    3) plot the results
    """
    num_plot_radius=30
    num_plot_speed=30


    smpling1 = np.linspace(-0.5, 0.5, num_plot_radius)
    smpling2 = np.linspace(-0.5, 0.5, num_plot_speed)

    samples_radius=range_radius*smpling1+mean_radius
    samples_speed=range_speed*smpling2+mean_speed

    h_plot=test_coef(coef, smpling1, smpling2, dis1, dis2, alf, beta)
    p_plot=test_coef(coef_p, smpling1, smpling2, dis1, dis2, alf, beta)
    samples_radius, samples_speed = np.meshgrid(samples_radius, samples_speed)

    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(samples_radius, samples_speed, h_plot, rstride=1, cstride=1, cmap='hot')

    ax.set_xlabel(name[0])
    ax.set_ylabel(name[1])
    ax.set_zlabel('height')
    plt.title('Height as a function of two input random variables')



    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(samples_radius, samples_speed, p_plot, rstride=1, cstride=1, cmap='hot')

    ax.set_xlabel(name[0])
    ax.set_ylabel(name[1])
    ax.set_zlabel('Particle mass flux')
    plt.title('Particle mass flux as a function of two input random variables')


    """
    Then resample random space with a huge amout of sample points (tensor product of sample points of two independent sample pooints in each direction)
    plug sample point into the expression
    got distribution of h
    """
    number_of_bins=100
    number_of_smaples_each_rv=500
    
    s1=generate_samples(dis1, number_of_smaples_each_rv, alf, beta)
    s2=generate_samples(dis2, number_of_smaples_each_rv, alf, beta)

    fig = plt.figure()
    h_hist=test_coef(coef, s1, s2, dis1, dis2, alf, beta)
    n, bins, patches = plt.hist(h_hist.flatten(), number_of_bins, normed=1, facecolor='green')
    plt.title('Histogram for Eruption Height')
    plt.show()

    fig = plt.figure()
    p_hist=test_coef(coef_p, s1, s2, dis1, dis2, alf, beta)
    n, bins, patches = plt.hist(p_hist.flatten(), number_of_bins, normed=1, facecolor='blue')
    plt.title('Histogram for Particle Mass Flux')
    plt.show()

    return 0
    
alf=0.5
beta=0.5
num_sample=[2,2]
mean_rv=[200, 200]
range_rv=[100, 100]
rv_value=[222, 333]
v_dis=["Gamma", "Gamma"]
v_name=['radius', 'speed']
PCQ_UQ(num_sample, mean_rv, range_rv, v_dis, v_name, rv_value, alf, beta)