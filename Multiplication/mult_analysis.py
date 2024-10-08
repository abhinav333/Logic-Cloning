# -*- coding: utf-8 -*-
"""
Analysis of error of approximate multiplication circuits. Computation of TP values of BW and Booth circuits.
Author: Abhinav
"""

# ctypes_test.py
import ctypes
import pathlib
import os
import numpy as np
import matplotlib.pyplot as mp

from numpy.ctypeslib import ndpointer


#Configure the bit-width
BIT_WIDTH=8

UNSIGNED_LOWER_LIMIT=0
UNSIGNED_UPPER_LIMIT=(1<<BIT_WIDTH)-1

# Range of signed number for bit width B is from -(2^(B-1)) to (2^(B-1))-1
SIGNED_HALF=1<<(BIT_WIDTH-1)
SIGNED_UPPER_LIMIT=((SIGNED_HALF)-1)
SIGNED_LOWER_LIMIT=-1*((SIGNED_HALF))


PP_NUM_SIGNED=((BIT_WIDTH//2))




def unsigned_to_signed(a,bit_width):
    b=0
    mask=pow(2,bit_width)
    if((a >> (bit_width-1) & 1)):
        b=-1*abs(a-mask)
    else:
        b=a
    return b

ARRAY_SIZE=PP_NUM_SIGNED*(BIT_WIDTH+2)

#ctypes initialization
if __name__ == "__main__":
    # Load the shared library into ctypes
    libname = str(os.path.join(pathlib.Path().absolute(),"libapp_8.so"))
    c_lib = ctypes.CDLL(libname)

#ctypes interface with C multiplication circuit implementations
array_1d_c=ctypes.POINTER(ctypes.c_int8)
array_1d_uint8= np.ctypeslib.ndpointer(dtype=np.int8, ndim=2, flags='CONTIGUOUS')

c_lib.unsigned_to_signed.argtypes = (ctypes.c_longlong,ctypes.c_ubyte)
c_lib.unsigned_to_signed.restype =  (ctypes.c_longlong)

c_lib.baugh_wooley_hdl.argtypes = (ctypes.c_longlong,ctypes.c_longlong)
c_lib.baugh_wooley_hdl.restype =  (ctypes.c_longlong)

c_lib.area_optimized_hdl_opt.argtypes = (ctypes.c_int64,ctypes.c_int64)
c_lib.area_optimized_hdl_opt.restype =  ctypes.c_int64

c_lib.area_optimized_hdl_approx.argtypes = (ctypes.c_int64,ctypes.c_int64)
c_lib.area_optimized_hdl_approx.restype = (ctypes.c_int64)


c_lib.lc_baugh_wooley_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_baugh_wooley_opt.restype = (ctypes.c_int64)

c_lib.lc_booth_opt.argtypes = (ctypes.c_int64,ctypes.c_int64,ctypes.c_int)
c_lib.lc_booth_opt.restype = (ctypes.c_int64)


c_lib.probe_values.argtypes = (array_1d_uint8,)
c_lib.probe_values.restype =  None

c_lib.probe_lut_result.argtypes = (array_1d_uint8,)
c_lib.probe_lut_result.restype =  None


array=np.zeros((PP_NUM_SIGNED,(BIT_WIDTH+2)),dtype=np.int8)
array_sum=np.zeros((PP_NUM_SIGNED,(BIT_WIDTH+2)))


array_val=np.zeros((BIT_WIDTH,(2*BIT_WIDTH)),dtype=np.int8)
array_val_sum=np.zeros((BIT_WIDTH,(2*BIT_WIDTH)))

error_prob_accurate=np.zeros(2*BIT_WIDTH)
error_prob_approx=np.zeros(2*BIT_WIDTH)

#Computation of TP values
def compute_probabilities(num,n_bits):
    error_prob=np.array(list(np.binary_repr(num,width=n_bits)),dtype=np.float64)
    return error_prob
    


max_error=0
max_errors=[]
max_error_multiplicand=[]
max_error_mutiplier=[]

range_mul_acc=[]  
range_mul_app=[]  

lut1_muc=[]
lut2_muc=[]
lut3_muc=[]
lut4_muc=[]
errors_series=[]
range_val=0

for app_mulc in range(SIGNED_LOWER_LIMIT, SIGNED_UPPER_LIMIT+1):
    for app_mult in range(SIGNED_LOWER_LIMIT, SIGNED_UPPER_LIMIT+1):
        range_val=range_val+1
        res1=app_mulc * app_mult
        range_mul_acc.append(res1)
        #res2=c_lib.baugh_wooley_hdl(ctypes.c_longlong(app_mulc),ctypes.c_longlong(app_mult))  
        #res2=c_lib.area_optimized_hdl_opt(ctypes.c_longlong(app_mulc),ctypes.c_longlong(app_mult))  
        res2=c_lib.lc_baugh_wooley_opt(ctypes.c_longlong(app_mulc),ctypes.c_longlong(app_mult),ctypes.c_int(2))
        #res2=c_lib.lc_booth_opt(ctypes.c_longlong(app_mulc),ctypes.c_longlong(app_mult),ctypes.c_int(1))
        #res2=c_lib.area_optimized_hdl_approx(ctypes.c_longlong(app_mulc),ctypes.c_longlong(app_mult))  
        c_lib.probe_lut_result(array)
        c_lib.probe_values(array_val)
        array_val_sum=array_val_sum+array_val
        array_sum=array_sum+array
        range_mul_app.append(res2)
        error_sign=np.int64(res1)-np.int64(res2)
        error=np.abs(error_sign)
        errors_series.append(error_sign)
        error_prob_accurate+=compute_probabilities(res1,2*BIT_WIDTH)
        error_prob_approx+=compute_probabilities(res2,2*BIT_WIDTH)
        if(max_error < error):
            max_error=error
            max_errors.append(error)
            max_error_multiplicand.append(app_mulc)
            max_error_mutiplier.append(app_mult)



error_prob_accurate=error_prob_accurate/range_val
error_prob_approx=error_prob_approx/range_val

array_val_sum=np.fliplr(array_val_sum)
array_sum=np.fliplr(array_sum)

array_booley_prob=array_val_sum/range_val
array_booth_prob=array_sum/range_val


