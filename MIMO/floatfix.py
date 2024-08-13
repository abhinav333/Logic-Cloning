#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 10:07:49 2021

@author: Abhinav
"""
import numpy as np
import ctypes
import pathlib
import os
from numpy.ctypeslib import ndpointer

libname = str(os.path.join(pathlib.Path().absolute(),"libapp.so"))
c_lib = ctypes.CDLL(libname)


BOOTH_APPROX=0
MUL8_14_4=1
MUL8_13_4=2
LC_BW_2=3
LC_BW_1=4
LC_BOOTH_2=5
LC_BOOTH_1=6
AXBM1=7
AXBM2=8



c_lib.area_optimized_hdl_approx.argtypes = (ctypes.c_int64,ctypes.c_int64)
c_lib.area_optimized_hdl_approx.restype = (ctypes.c_int64)

c_lib.multilevel_signed_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_char,ctypes.c_char)
c_lib.multilevel_signed_opt.restype = (ctypes.c_int64)

c_lib.multilevel2B_signed_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_char,ctypes.c_char)
c_lib.multilevel2B_signed_opt.restype = (ctypes.c_int64)

c_lib.lc_baugh_wooley_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_baugh_wooley_opt.restype = (ctypes.c_int64)

c_lib.lc_baugh_wooley_opt_ac.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_baugh_wooley_opt_ac.restype = (ctypes.c_int64)

c_lib.lc_booth_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_booth_opt.restype = (ctypes.c_int64)

c_lib.lc_booth_opt_ac.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_booth_opt_ac.restype = (ctypes.c_int64)

c_lib.ax_bm1_opt.argtypes = (ctypes.c_int64,ctypes.c_int64)
c_lib.ax_bm1_opt.restype = (ctypes.c_int64)

c_lib.ax_bm2_opt.argtypes = (ctypes.c_int64,ctypes.c_int64)
c_lib.ax_bm2_opt.restype = (ctypes.c_int64)



    

def approximate_multiplier(multiplicand, multiplier,m_config):
    if(m_config==BOOTH_APPROX):
        result=c_lib.area_optimized_hdl_approx(ctypes.c_int64(multiplicand),ctypes.c_int64(multiplier));
    elif(m_config==MUL8_13_4):
        result=result=c_lib.multilevel_signed_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_char(3),ctypes.c_char(4))
    elif(m_config==MUL8_14_4):
        result=result=c_lib.multilevel_signed_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_char(4),ctypes.c_char(4))
    elif(m_config==LC_BW_2):
        result=c_lib.lc_baugh_wooley_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_int(2));
    elif(m_config==LC_BW_1):
        result=c_lib.lc_baugh_wooley_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_int(1));
    elif(m_config==LC_BOOTH_2):
        result=c_lib.lc_booth_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_int(2));
    elif(m_config==LC_BOOTH_1):
        result=c_lib.lc_booth_opt(ctypes.c_int64(multiplicand), ctypes.c_int64(multiplier),ctypes.c_int(1));
    elif(m_config==AXBM1):
        result=c_lib.ax_bm1_opt(ctypes.c_int64(multiplicand),ctypes.c_int64(multiplier));
    elif(m_config==AXBM2):
        result=c_lib.ax_bm2_opt(ctypes.c_int64(multiplicand),ctypes.c_int64(multiplier));
    else:
        print("Invalid choice\n");
    return(result)    



    
def float2fix_scale(a,m,n):  # m is the bit-width, n is fractional part
    scaled=np.multiply(a,np.power(2,n))
    if(scaled==np.NaN):
        print("Nan detected")
    scaled=np.int(round(scaled))
    min_bound=-1*np.power(2,m-1)
    max_bound=np.power(2,m-1)-1/np.power(2,1*n)
    return np.int(np.clip(scaled,min_bound,max_bound))

def float2fix_double_descale(a,m,n):  # m is the bit-width, n is fractional part
    descaled=np.divide(a,np.power(2,2*n))
    min_bound=-1*np.power(2,m-1)
    max_bound=np.power(2,m-1)-1/np.power(2,1*n)
    return np.clip(descaled,min_bound,max_bound)

def float2fix_descale(a,m,n):  # m is the bit-width, n is fractional part
    descaled=np.divide(a,np.power(2,n))
    min_bound=-1*np.power(2,m-1)
    max_bound=np.power(2,m-1)-1/np.power(2,1*n)
    return np.clip(descaled,min_bound,max_bound)



##Change the multiplier configuration here    
def sfix_mul(a,b,m,n,multiplier_config,bit_width):
    op_a=float2fix_scale(a, m, n)
    op_b=float2fix_scale(b, m, n)
    op_c=approximate_multiplier(op_a, op_b ,multiplier_config)
    c=float2fix_double_descale(op_c, m, n)
    return c



def sfix_mul_cmplx(a,b,m,n,multiplier_config):
    c_real_1=sfix_mul(a.real,b.real,m,n,multiplier_config)
    c_real_2=sfix_mul(a.imag,b.imag,m,n,multiplier_config)
    c_real=c_real_1-c_real_2
    c_imag_1=sfix_mul(a.imag,b.real,m,n,multiplier_config)
    c_imag_2=sfix_mul(a.real,b.imag,m,n,multiplier_config)
    c_imag=c_imag_1+c_imag_2
    return np.complex128(c_real+c_imag*1j)
    
def sfix_div_cmplx(a,b,m,n,multiplier_config):
    b_conj=np.conjugate(b)
    num_1=sfix_mul_cmplx(a,b_conj,m,n,multiplier_config)
    b_real_2=sfix_mul(b.real,b.real,m,n,multiplier_config)
    b_imag_2=sfix_mul(b.imag,b.imag,m,n,multiplier_config)
    b_2=b_real_2+b_imag_2
    c_real=num_1.real/b_2
    c_imag=num_1.imag/b_2
    return np.complex128(c_real+c_imag*1j)



def sfix_mat_mul_cmplx(A,B,m,n,multiplier_config):
    if(A.shape[1]!=B.shape[0]):
        print("Matrix multiply error:dimension mismatch")
        return -1;
    C=np.zeros((A.shape[0],B.shape[1]), dtype=A.dtype)
    for row in range(0, A.shape[0]):
        for col in range(0,B.shape[1]):
            for su in range(0,A.shape[1]):
                p=sfix_mul_cmplx(A[row][su],B[su][col],m,n,multiplier_config)
                C[row][col]=C[row][col]+p
    return C             
    
    
#
def sfix_mat_inverse(A_mat,m,n,multiplier_config):
#A=np.random.randn(3,3)
    p=A_mat.shape[0]
    AI=np.eye(p,dtype=A_mat.dtype)
    A=np.block([A_mat,AI])
   
    # Applying Guass Jordan Elimination
    for i in range(p):
        if np.abs(A[i][i]) == 0.0:
            print('Divide by zero detected!')
        
        for j in range(p):
            if i != j:
                #ratio = A[j][i]/A[i][i]
                ratio = sfix_div_cmplx(A[j][i], A[i][i], m, n, multiplier_config)    
                for k in range(2*p):
                    factor=sfix_mul_cmplx(ratio, A[i][k], m, n, multiplier_config)
                    A[j][k] = A[j][k] - factor
    
    # Row operation to make principal diagonal element to 1
    for i in range(p):
        divisor = A[i][i]
        for j in range(2*p):
            A[i][j] = sfix_div_cmplx(A[i][j], divisor, m, n, multiplier_config)
    
    return A[:,p:]        
    



