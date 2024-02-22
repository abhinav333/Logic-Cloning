#ifndef _APP_
#define _APP_

#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define BIT_MSB1 3
#define BIT_MSB 2   //2i+1
#define BIT_CSB 1   //2i
#define BIT_LSB 0   //2i-1

#define BIT_WIDTH 16
#define PP_NUM_UNSIGNED ((BIT_WIDTH/2)+1)
#define PP_NUM_SIGNED ((BIT_WIDTH/2))
#define NEGATIVE_ZERO (1<<BIT_WIDTH+1)-1

#define LUT_TEST 99
#define LUT_AO_A 0
#define LUT_AO_B 1
#define LUT_AO_Z 2
#define LUT_AO_A1 3
#define LUT_AO_A2 4
#define LUT_AO_Am 5

#define LUT_ES_A 6
#define LUT_ES_B 7
#define LUT_ES_C 8

#define LUT_AXBM_A 9
#define LUT_AXBM_EN_1 10
#define LUT_AXBM_EN_2 11
#define LUT_BW_A 12
#define LUT_BW_B 13
#define LUT_BW_C 14



#define LUT_BWO_J 15
#define LUT_BWO_K 16
#define LUT_BWO_L 17
#define LUT_BWO_M 18
#define LUT_BWO_N 19


#define NMED 2
#define MRED 1
#define MED 0




#define UNSIGNED_LOWER_LIMIT 0
#define UNSIGNED_UPPER_LIMIT (1<<BIT_WIDTH)-1

// Range of signed number for bit width B is from -(2^(B-1)) to (2^(B-1))-1
#define SIGNED_HALF 1<<(BIT_WIDTH-1)
#define SIGNED_UPPER_LIMIT ((SIGNED_HALF)-1)
#define SIGNED_LOWER_LIMIT -1*((SIGNED_HALF))


#define BIT_WIDTH_2B 2*BIT_WIDTH
#define UNSIGNED_LOWER_LIMIT_2B 0
#define UNSIGNED_UPPER_LIMIT_2B (1<<BIT_WIDTH_2B)-1

// Range of signed number for bit width B is from -(2^(B-1)) to (2^(B-1))-1
#define SIGNED_HALF_2B 1<<(BIT_WIDTH_2B-1)
#define SIGNED_UPPER_LIMIT_2B ((SIGNED_HALF_2B)-1)
#define SIGNED_LOWER_LIMIT_2B -1*((SIGNED_HALF_2B))




#define N_RANGE (SIGNED_UPPER_LIMIT-SIGNED_LOWER_LIMIT+1)^2




#define CARRY_PROBE


#if BIT_WIDTH==8
#define PP_NUM_UNSIGNED_AXBM 3 //for 8-bit operand.
#elif BIT_WIDTH==16
#define PP_NUM_UNSIGNED_AXBM 6 //for 16-bit operand.
#elif BIT_WIDTH==17
#define PP_NUM_UNSIGNED_AXBM 6 //for 16-bit operand.
#elif BIT_WIDTH==4
#define PP_NUM_UNSIGNED_AXBM 2 //for 16-bit operand.
#elif BIT_WIDTH==5
#define PP_NUM_UNSIGNED_AXBM 1 //for 16-bit operand.
#endif

#define BIT_WIDTH_ML 8
#define BIT_WIDTH_2ML 16

//MULTIPLIERS
#define BOOTH_APPROX 0
#define MUL8_14_4 1
#define MUL8_13_4 2
#define LC_BW_2 3
#define LC_BW_2_AC 4
#define LC_BOOTH_2 5
#define LC_BOOTH_2_AC 6
#define AXBM1 7
#define AXBM2 8




typedef unsigned char uchar;
void print_binary(unsigned long long n,unsigned int k_bits);

void print_binary2(unsigned long long n,unsigned int k_bits);



void CARRY1(uchar s,uchar d, uchar cin, uchar* cout, uchar* xor_cout);


//Cin has to be divided into CYINIT & CI ::Note this is for VHDL implementation. You can do this for VHDL impelemntation.
void CARRY4(uchar* s, uchar* d, uchar Cin, uchar* cout_a, uchar* xor_cout_a);

uchar not_gate(uchar number);

//Function to create 64 bit value from a boolean function
//Area-Optimized Accurate and Approximate Softcore Signed Multiplier Architectures : CODE : AO
//An Efficient Softcore Multiplier Architecture for Xilinx FPGAs: CODE : ES

void LUT_6_2(uchar a, uchar b, uchar c, uchar d, uchar e, uchar f, uchar lut_type, int* counter, uchar* result);

long long approximate_multiplier(long long multiplicand, long long multiplier, uchar MULT);
void lut_testing();



/********************************
            HA - HA - HA - HA
            FA - FA - FA
Generate this kind of structure.
Now create strings of the oprands to be passed to each of these rows.

Row1: {op1, op2,.....}
Row2: {op1,op2.......}

like this where op1, op2... belong to Partial Product bits.
********************************/
/*****************************************
An Efficient Softcore Multiplier Architecture for Xilinx FPGAs
Martin Kumm, Shahid Abbas and Peter Zipf
University of Kassel, Germany
Digital Technology Group
******************************************/

long efficient_softcore(long multiplicand, long multiplier);



long efficient_softcore_hdl(long multiplicand, long multiplier);



long long unsigned_to_signed(unsigned long long a,uchar bit_width);

//signed

long long area_optimized_hdl(long long multiplicand, long long multiplier);

long long area_optimized_hdl_opt(long long multiplicand, long long multiplier);
long long area_optimized_hdl_opt_v2(long long multiplicand, long long multiplier);
long long area_optimized_hdl_opt_v3(long long multiplicand, long long multiplier);

void probe_carry(uchar *carry_values);
void probe_lut_result(uchar *lut_values);



//multplier definations
long long lc_baugh_wooley_opt(long long multiplicand, long long multiplier, int version);
long long lc_baugh_wooley_opt_ac(long long multiplicand, long long multiplier, int version);
long long area_optimized_hdl_approx(long long multiplicand, long long multiplier);
long long lc_booth_opt(long long multiplicand, long long multiplier,int version);
long long lc_booth_opt_ac(long long multiplicand, long long multiplier,int version);
long long multilevel_signed_opt(long long multiplicand_sign, long long multiplier_sign, char c2, char k_paper);
long long ax_bm1_opt(long multiplicand, long multiplier);
long long ax_bm2_opt(long multiplicand, long multiplier);


long long signed_multiplication(long long multiplicand, long long multiplier);


/*************
AxBM
************/




long long ax_bm1_hdl(long multiplicand, long multiplier);



void error_benchmarks(double* errors);





void error_compute();

void error_compute_signed();




#endif
