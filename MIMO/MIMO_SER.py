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

symerr_all=np.zeros((5,11,snr_iter.size))

scale=np.sqrt(cluster_number)

for K in range(K,K+2,2):
    H=np.random.normal(loc=0, scale=1, size=(M,K*2)).view(np.complex128)/np.sqrt(2*M)
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


pl.title('Symbol error rate: U={} B={} C={} QAM={} Iter={}'.format(K,M,cluster_number,M_QAM,set_iteration)) 
pl.xlabel('SNR')
pl.ylabel('Symbol error rate')
pl.legend()



pl.grid(True)



