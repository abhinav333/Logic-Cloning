 # -*- coding: utf-8 -*-6464
"""
Created on Thu May 16 11:29:45 2019

@author: Abhinav_Kulkarni
MIMO uplink decoding codes
"""
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as pl

import floatfix_8 as fp        


BOOTH_APPROX=0
MUL8_14_4=1
MUL8_13_4=2
LC_BW_2=3
LC_BW_2_AC=4
LC_BOOTH_2=5
LC_BOOTH_2_AC=6
AXBM1=7
AXBM2=8
'''
Notes of simulation patters:
    
    1. 8 bit  -> shift by 7, 4 qam 



of it, either on the Editor or the Console.






'''

#channel_matrix=sp.loadmat("H_matrix_128.mat")

#H_scm=channel_matrix['H_mat']


#function declarations
#QPSK Modulation
def qamod(sym,constel):  
    c=np.sqrt(constel)
    re=-2*np.mod(sym,c)+c-1
    im=2*np.floor(sym/c)-c+1
    return re+1j*im

def qademod(symb,constel):  #One symbol at a time
    QAM_map=qamod(np.arange(0,constel),constel)
    symbol_distance=np.square(np.abs(symb-QAM_map))
    return np.argmin(symbol_distance)

def calculate_evm(ex_est,ex):     
    evm_real=np.sum(np.square(ex_est.real-ex.real))
    evm_imag=np.sum(np.square(ex_est.imag-ex.imag))
    return np.sqrt((evm_real+evm_imag)/ex.size)
        
def calculate_symbol_error_rate(x_est,data): #requires a flattened array
    error_rate=0
    for loop in range(data.size):
        if qademod(x_est[loop],M_QAM)!=data[loop]:
            error_rate+=1
    return (error_rate)
        
def calculate_SINR_AVG(x_est_l,x_l):  #do not put x_est_frame = x_frame, else division by zero error 
    power_transmit=np.square(np.linalg.norm(x_l,axis=0))
    power_receive=np.square(np.linalg.norm(x_est_l,axis=0))
    return np.average(10*np.log10(np.divide(power_receive,np.abs((power_transmit-power_receive)))))
#    return np.average(np.abs(power_transmit/power_transmit-power_receive))

def calculate_SigmaN2(SNR):
    return np.float_power(10,-1*SNR/10.0) #evaluates to float64



#variable declarations
#MKratio=32   #MK ratio
set_iteration=5
K=np.uint16(4)    #number of user t  erminals
cluster_number=np.uint16(4) #number of clusters  
MperC=32
M=np.uint16(cluster_number*MperC)   #number of antennas 

bit_width=8
frac_width=7
multiplier_config=1

multiple_frame=5
M_QAM=np.uint16(4)
Nframes=1

mmse_roh=0.01

QAM_map=qamod(np.arange(0,M_QAM),M_QAM)     #write a qammod function
data= np.random.randint(0,M_QAM,size=(K,Nframes))
QAM_mod=np.array(QAM_map[data]) #symbol mapped tx signal
QAM_var=np.sqrt(np.var(QAM_map))
x=QAM_mod
x_d=x

def generate_data(M_QAM,K,Nframes):
    QAM_map= qamod(np.arange(0,M_QAM),M_QAM)     #write a qammod function
    data= np.random.randint(0,M_QAM,size=(K,Nframes))
    QAM_mod=np.array(QAM_map[data]) #symbol mapped tx signal
    QAM_var=np.sqrt(np.var(QAM_map))
    x=QAM_mod/QAM_var
    #x=QAM_mod
    return x,data,QAM_var

def simulate_receive_signal(x,SNR,scale,H):
    var=np.sqrt(calculate_SigmaN2(SNR)/2)
    noise_sys=var*np.random.normal(loc=0, scale=1, size=(M,Nframes*2)).view(np.complex128)
    y=np.zeros((M,Nframes),dtype=x.dtype) #array holding y(number of antennas)
    for f in range(0,Nframes):
        y[:,f]= (H @ x[:,f]) + noise_sys[:,f]
    return y        
     
x_est=np.zeros(x.shape,dtype=x.dtype)  
########################################################   
##ADMM-gauss siedel

i_gauss_ADMM=np.uint16(6) #number of iterations

#ZF method

i_gauss_ADMM=np.uint16(6) #number of iterations

#ZF method
def zf_method(y,x,A):
    x_est_zf=np.zeros(x.shape,dtype=x.dtype) 
    invA=np.zeros(A.shape,dtype=A.dtype)
    yMF=np.zeros((K,Nframes),dtype=x.dtype)
    invA=np.linalg.inv(A)    
    for f in range(0,Nframes):
        yMF[:,f]=np.conjugate(H.T) @ y[:,f]   #read about the MRC ratio combining
        x_est_zf[:,f]=invA @ yMF[:,f]
    return x_est_zf

def zf_method_mat(y,x,A,multiplier_config):
    x_est_zf=np.zeros(x.shape,dtype=x.dtype) 
    invA=np.zeros(A.shape,dtype=A.dtype)
    yMF=np.zeros((K,Nframes),dtype=x.dtype)
    invA=fp.sfix_mat_inverse(A, bit_width, frac_width, multiplier_config)    #
    for f in range(0,Nframes):
        yMF[:,[f]]=fp.sfix_mat_mul_cmplx(np.conjugate(H.T), y[:,[f]], bit_width, frac_width, multiplier_config)   #read about the MRC ratio combining
        x_est_zf[:,[f]]=fp.sfix_mat_mul_cmplx(invA, yMF[:,[f]], bit_width, frac_width, multiplier_config) 
    return x_est_zf



def mmse_method(y,x,A,sigmaN2):
    x_est_mmse=np.zeros(x.shape,dtype=x.dtype) 
    invA=np.zeros(A.shape,dtype=A.dtype)
    yMF=np.zeros((K,Nframes),dtype=x.dtype)
    sigmaN2_QAM=np.square(QAM_var)
    A_mmse=A+((sigmaN2/sigmaN2_QAM)*np.identity(K,dtype=np.complex128))
    invA=np.linalg.inv(A_mmse)    
    for f in range(0,Nframes):
        yMF[:,f]=np.conjugate(H.T) @ y[:,f]   #read about the MRC ratio combining
        x_est_mmse[:,f]=invA @ yMF[:,f]
    return x_est_mmse

snr_iter=np.arange(-5,30,2)

#time_admm_r1=np.zeros(snr_iter.size)
#time_admm_r2=np.zeros(snr_iter.size)
#time_admm_r3=np.zeros(snr_iter.size)
#symerr_admm_r1=np.zeros(snr_iter.size) 
#symerr_admm_r2=np.zeros(snr_iter.size) 

#symerr_admm_r3=np.zeros(snr_iter.size)  
#time_cd=np.zeros(snr_iter.size)
#symerr_cd=np.zeros(snr_iter.size) 
#time_sgd=np.zeros(snr_iter.size)
#symerr_sgd=np.zeros(snr_iter.size) 
#time_zf=np.zeros(snr_iter.size)
symerr_zf=np.zeros(snr_iter.size) 
symerr_zf_ba=np.zeros(snr_iter.size) 
symerr_zf_m13_2=np.zeros(snr_iter.size) 
symerr_zf_m14_2=np.zeros(snr_iter.size) 
symerr_zf_lc_bw_2=np.zeros(snr_iter.size) 
symerr_zf_lc_bw_ac_2=np.zeros(snr_iter.size) 
symerr_zf_lc_booth_2=np.zeros(snr_iter.size) 
symerr_zf_lc_booth_ac_2=np.zeros(snr_iter.size) 
symerr_zf_axbm1=np.zeros(snr_iter.size)
symerr_zf_axbm2=np.zeros(snr_iter.size)
#symerr_bsgd=np.zeros(snr_iter.size) 
#symerr_bsgd1=np.zeros(snr_iter.size) 
#symerr_bsgd2=np.zeros(snr_iter.size)
#symerr_amp_fd=np.zeros(snr_iter.size) 
#symerr_amp_pd=np.zeros(snr_iter.size)
#symerr_ep=np.zeros(snr_iter.size)  
#symerr_mmse=np.zeros(snr_iter.size) 

symerr_all=np.zeros((5,11,snr_iter.size))

scale=np.sqrt(cluster_number)

for K in range(K,K+2,2):
    #print('U={}'.format(K))
    H=np.random.normal(loc=0, scale=1, size=(M,K*2)).view(np.complex128)/np.sqrt(2*M)
   # H=H_scm/np.sqrt(M)
    #H=2*H
    A=np.transpose(np.conjugate(H)) @ H   #zero forcing
    mmse_roh=0.01
    beta=np.float64(K/M)
    A_mmse=A+mmse_roh*np.identity(K,dtype=np.complex128)
    d_size=K*Nframes*multiple_frame
    for i in range(snr_iter.size):
       # print(snr_iter[i])
        for fr in range(multiple_frame):
           # print(fr)
            print('U={} SNR={} FRAME={}'.format(K,snr_iter[i],fr))
            x=generate_data(M_QAM,K,Nframes)  
            y=simulate_receive_signal(x[0],snr_iter[i],scale,H)
            x_est=zf_method(y,x[0],A)
            symerr_zf[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,BOOTH_APPROX)
            symerr_zf_ba[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,MUL8_13_4)
            symerr_zf_m13_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,MUL8_14_4)
            symerr_zf_m14_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,LC_BW_2)
            symerr_zf_lc_bw_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,LC_BW_2_AC)
            symerr_zf_lc_bw_ac_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,LC_BOOTH_2)
            symerr_zf_lc_booth_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,LC_BOOTH_2_AC)
            symerr_zf_lc_booth_ac_2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,AXBM1)
            symerr_zf_axbm1[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
            x_est=zf_method_mat(y,x[0],A,AXBM2)
            symerr_zf_axbm2[i]+=calculate_symbol_error_rate(x[2]*x_est.flatten(),x[1].flatten())
        symerr_zf[i]=symerr_zf[i]/d_size 
        symerr_zf_ba[i]=symerr_zf_ba[i]/d_size 
        symerr_zf_m13_2[i]=symerr_zf_m13_2[i]/d_size 
        symerr_zf_m14_2[i]=symerr_zf_m14_2[i]/d_size 
        symerr_zf_lc_bw_2[i]=symerr_zf_lc_bw_2[i]/d_size 
        symerr_zf_lc_bw_ac_2[i]=symerr_zf_lc_bw_ac_2[i]/d_size 
        symerr_zf_lc_booth_2[i]=symerr_zf_lc_booth_2[i]/d_size
        symerr_zf_lc_booth_ac_2[i]=symerr_zf_lc_booth_ac_2[i]/d_size
        symerr_zf_axbm1[i]=symerr_zf_axbm1[i]/d_size
        symerr_zf_axbm2[i]=symerr_zf_axbm2[i]/d_size
        
#        symerr_bsgd[i]=symerr_bsgd[i]/d_size 
#        symerr_bsgd1[i]=symerr_bsgd1[i]/d_size 
#        symerr_amp_fd[i]=symerr_amp_fd[i]/d_size 
#        symerr_amp_pd[i]=symerr_amp_pd[i]/d_size
#        symerr_ep[i]=symerr_ep[i]/d_size  
#        symerr_mmse[i]=symerr_mmse[i]/d_size
    symerr_all[K//2,:,:]=np.array([snr_iter,symerr_zf,symerr_zf_ba,symerr_zf_m13_2,symerr_zf_m14_2,symerr_zf_lc_bw_2,symerr_zf_lc_bw_ac_2,symerr_zf_lc_booth_2,symerr_zf_lc_booth_ac_2,symerr_zf_axbm1,symerr_zf_axbm2])
     
     
    



pl.semilogy(snr_iter,symerr_zf_axbm2,'r-x',label='AXBM2'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_axbm1,'r-x',label='AXBM1'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_lc_booth_ac_2,'b-s',label='LC BOOTH AC'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_lc_booth_2,'m-^',label='LC BOOTH')
pl.semilogy(snr_iter,symerr_zf_lc_bw_ac_2,'g-o',label='LC BW AC'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_lc_bw_2,'y--p',label='LC BW'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_m14_2,'c--p',label='MUL8_14_2'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_m13_2,'c--^',label='MUL8_13_2'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf_ba,'y--o',label='BOOTH APPROX'.format(set_iteration))
pl.semilogy(snr_iter,symerr_zf,'k-*',label='ZF',alpha=0.7)
#pl.semilogy(snr_iter,symerr_mmse,'y-s',label='MMSE',alpha=0.2)


pl.title('Symbol error rate: U={} B={} C={} QAM={} Iter={}'.format(K,M,cluster_number,M_QAM,set_iteration)) 
pl.xlabel('SNR')
pl.ylabel('Symbol error rate')
pl.legend()


""" 
pl.plot(QAM_var*(x_est.flatten().real),QAM_var*(x_est.flatten().imag),'bo')
pl.plot(QAM_map.real,QAM_map.imag,'ro')
pl.title('Co-ordinate descent: SNR={}dB u-parameter={}')  
pl.xlabel('Real Axis')
pl.ylabel('Imaginary Axis')
"""


pl.grid(True)







#pl.plot(admm_i_gauss_ADMM,evm_simulation,'b--o')  
#pl.title('ADMM-GS:  SNR={}dB   User-Terminals={}    Antenna/UT={}    Clusters={}'.format(SNR,K,MKratio,Max_C))  
#pl.xlabel('Number of iterations')
#pl.ylabel('Error Vector Magnitude (EVM)')
#pl.grid(True)
#print (calculate_evm(x_est_admm,x))

##################################################################################################

   
#Extra stuff
       
##ZF method
#invA=np.zeros(A.shape,dtype=A.dtype)
#invA=np.linalg.inv(A)    
#for f in range(0,Nframes):
#invA=np.zeros(A.shape,dtype=A.dtype)
#invA=np.linalg.inv(A)    
#    x_est[:,f]=invA @ yMF[:,f]
#
#x_est_inv=np.copy(x_est)                


#Gauss siedel
#estimate the initial guess from by x=inverse(A)b. As A is diagonally dominant inverse of A can be obtained by reciprocal of the A. 
#Aid=np.zeros(A.shape,dtype=A.dtype)
#np.fill_diagonal(Aid,(1/np.diagonal(A)))
#n_iter=3
 #number of iteration to perform 
#try for 1st frame  
#x_init=Aid @ yMF[:,0]
#x_init_t=np.copy(x_init)
#for i in range (n_iter):
#    for r_scan in range (K):
#        x_init_t[r_scan]=yMF[r_scan,0]
#        for c_scan in range (K):
#            if r_scan!=c_scan:
#                x_init_t[r_scan]-=(A[r_scan,c_scan] * x_init[c_scan])
#        #c_scan loop
#        x_init[r_scan]=(1/A[r_scan,r_scan])*(x_init_t[r_scan])

#err=np.abs(x_est.reshape(-1,1)-x.reshape(-1,1))
#pl.plot(err,'b')




