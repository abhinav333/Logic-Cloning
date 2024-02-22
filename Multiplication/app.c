
//
//Author :Abhinav


#include "app.h"
#include "math.h"

#include "stdio.h"
#include "stdlib.h"

#define M1
#define MUM 2*M1



//global variables
uchar carry_opt[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
uchar lut_opt[PP_NUM_SIGNED][BIT_WIDTH+2]={ };

uchar pp_probe[BIT_WIDTH][2*BIT_WIDTH]={ };

void write_csv(double *array, int range,char filename[])
{
    FILE *fpt;
    fpt = fopen(filename, "w+");
    for(int loop=0;loop<range;loop++)
    {
       fprintf(fpt,"%lf\n",array[loop]);
    }
    fclose(fpt);

}

void compute_probabilities(long long number, double bit_prob[2*BIT_WIDTH])
{
    for(int sc=0;sc<(2*BIT_WIDTH);sc++)
    {
        bit_prob[sc]+=((number >> sc) & 1) ? 1 : 0;
    }
}
void noramlize_probabilities(double range,double bit_prob[2*BIT_WIDTH])
{
    for(int sc=0;sc<(2*BIT_WIDTH);sc++)
    {
        bit_prob[sc]/=range;
    }

}
void log_prob(double ln_prob_0[2*BIT_WIDTH],double ln_prob_1[2*BIT_WIDTH],double bit_prob[2*BIT_WIDTH])
{
    for(int sc=0;sc<(2*BIT_WIDTH);sc++)
    {
        ln_prob_0[sc]=log(1.0-bit_prob[sc]);
        ln_prob_1[sc]=log(bit_prob[sc]);

    }
}

double compute_probability_num(long long number,double ln_prob_0[2*BIT_WIDTH],double ln_prob_1[2*BIT_WIDTH])
{
    double prob=0;
    for(int sc=0;sc<(2*BIT_WIDTH);sc++)
    {
        prob+=((number >> sc) & 1) ? ln_prob_1[sc] : ln_prob_0[sc];
    }
    return(prob);
}

void compute_probability_distribution(double prob[1<<(2*BIT_WIDTH)],double ln_prob_0[2*BIT_WIDTH],double ln_prob_1[2*BIT_WIDTH])
{
        for(long long b=-(1<<((2*BIT_WIDTH)-1)); b<(1<<((2*BIT_WIDTH)-1));b++)
        {
            prob[b+(1<<((2*BIT_WIDTH)-1))]=compute_probability_num(b,ln_prob_0,ln_prob_1);
        }
}


double m_oldM, m_newM, m_oldS, m_newS;
int m_n;
void Push(double x)
{
            m_n++;
            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (m_n == 1)
            {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
            }
            else
            {
                m_newM = m_oldM + (x - m_oldM)/m_n;
                m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

                // set up for next iteration
                m_oldM = m_newM;
                m_oldS = m_newS;
            }
}

double Mean()
{
   return (m_n > 0) ? m_newM : 0.0;
}

double Variance()
{
   return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
}

double StandardDeviation()
{
   return sqrt( Variance() );
}


typedef unsigned char uchar;
void print_binary(unsigned long long n,unsigned int k_bits)
{
	uchar bit=0;
	for (int i=0;i<k_bits;i++)
	{
		bit= n & ((1<<(k_bits-1))>>i) ? 1 : 0 ;
		printf("%d",bit);
	}
}

void print_binary2(unsigned long long n,unsigned int k_bits)
{
	uchar bit=0;
	for (int i=k_bits-1;i>=0;i--)
	{
		bit= n & (1<<i) ? 1 : 0 ;
		printf("%d",bit);
	}
}

long long unsigned_to_signed(unsigned long long a,uchar bit_width)
{
    long long b=0;
    unsigned long long mask=pow(2,bit_width);
    if((a >> (bit_width-1) & 1))
        b=-1*abs(a-mask);
    else
        b=a;
    return b;
}


void CARRY1(uchar s,uchar d, uchar cin, uchar* cout, uchar* xor_cout)
{
    switch(s)
        {
            case 0:
                *cout=d;
                break;
            case 1:
                *cout=cin;
                break;
            default:
                printf("Invalid s\n");
                break;
        }
      *xor_cout=(cin ^ s) & 1;
}

//Cin has to be divided into CYINIT & CI ::Note this is for VHDL implementation. You can do this for VHDL impelemntation.
void CARRY4(uchar* s, uchar* d, uchar Cin, uchar* cout_a, uchar* xor_cout_a)
{
    uchar c=Cin;
    for(int i=0;i<4; i++)
    {
        CARRY1(s[i],d[i],c,cout_a[i],xor_cout_a[i]);
        c=cout_a[i];
    }
}

uchar not_gate(uchar number)
{
    return (~number & 1);
}

//Function to create 64 bit value from a boolean function
//Area-Optimized Accurate and Approximate Softcore Signed Multiplier Architectures : CODE : AO
//An Efficient Softcore Multiplier Architecture for Xilinx FPGAs: CODE : ES

void LUT_6_2(uchar a, uchar b, uchar c, uchar d, uchar e, uchar f, uchar lut_type, int* counter, uchar* result)
{
    uchar ls=0,lc=0,lz=0,lo=0,lg=0,lo1=0,lo2=0,lo3=0,lse=0,llo2=0,li=0,lll=0;
    uchar la1=0,la2=0,la3=0,l4x=0,lor1=0,lxor1=0,lx1=0,lx2=0;
    uchar sig1=0,sig2=0;
    switch(lut_type)
        {
            case LUT_AO_A:
                ls=(not_gate(a) & b & c) | (a & not_gate(b) & not_gate(c));
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lz=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & c);
                lo1=(not_gate(ls) & e) | (ls & f);
                lo2=(not_gate(lc) & lo1) | (lc & not_gate(lo1));
                lo3=(not_gate(lz) & lo2);
                lo=lo3 ^ d;
                lg=d;
                break;

            case LUT_AO_B:
                lse=(a & b & c) | (not_gate(a) & not_gate(b) & not_gate(c)) | (not_gate(a) & b & not_gate(e)) | (a & not_gate(b) & e) | (not_gate(a) & not_gate(b) & c & not_gate(e)) | (a & b & not_gate(c) & e);
                lo=lse ^ d;
                lg=d;
                break;

            case LUT_AO_Z:
                lo=1;
                break;

            case LUT_AO_A1:
                ls=(not_gate(a) & b & c) | (a & not_gate(b) & not_gate(c));
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lz=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & c);
                lo1=(not_gate(ls) & e) | (ls & f);
                lo2=(not_gate(lc) & lo1) | (lc & not_gate(lo1));
                lo3=(not_gate(lz) & lo2);
                //outputs
                llo2=not_gate(lz) & lc;
                lll=(lo3 ^ d);
                lo=lll ^ llo2;
                lg=(lll & llo2) | (d & not_gate(lll)); //Implementation of a mux (My claim)
                //lg=(lll & llo2); //Given in the paper.
                break;

            case LUT_AO_A2:
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lz=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & c);
                lo1=e;
                lo2=(not_gate(lc) & lo1) | (lc & not_gate(lo1));
                lo3=(not_gate(lz) & lo2);
                lo=lo3 ^ not_gate(d);
                lg=not_gate(d);
                break;

            case LUT_AO_Am:
                ls=(not_gate(a) & b & c) | (a & not_gate(b) & not_gate(c));
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lz=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & c);
                lo1=(not_gate(ls) & e) | (ls & f);
                lo2=(not_gate(lc) & lo1) | (lc & not_gate(lo1));
                lo3=(not_gate(lz) & lo2);

                li=(not_gate(a) & b & c ) | (a & not_gate(b)) | (a & b & not_gate(c));
                //outputs
                lo=(lo3 ^ d) ^ (li);
                lg=0;
                break;

            case LUT_ES_A:
                ls=(not_gate(a) & b & c) | (a & not_gate(b) & not_gate(c));
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lz=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & c);
                lo1=(not_gate(ls) & e) | (ls & f);
                lo2=(not_gate(lc) & lo1) | (lc & not_gate(lo1));
                lo3=(not_gate(lz) & lo2);
                //outputs
                lo=lo3 ^ d;
                lg=d;
                break;

            case LUT_ES_B:
                lc=(a & not_gate(b) )| (a & b & not_gate(c));
                lo3=not_gate(lc);
                //outputs
                lo=lo3 ^ d;
                lg=d;
                break;

            case LUT_ES_C:
                //outputs
                lo=1;
                lg=d;
                break;

            case LUT_AXBM_A:
                //a= x1 b=x2 c=s d=m2 e=m1 f=m0
                l4x=not_gate(a ^ b);
                lor1=(d & a) | (e & b) | (f & l4x);
                lxor1=(lor1 ^ c);
                lo=lxor1;
                lg=0;
                break;

            case LUT_AXBM_EN_1:
                //Verified OK.
                //a=di+2 b=di+1 c=d d=di-1
                lx1=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & d) | (not_gate(a) & not_gate(b) & c & not_gate(d))|(a & b & c & not_gate(d));
                lx2=(not_gate(b) & c & d) | (not_gate(a) & b & not_gate(c) & not_gate(d)) | (a & b & not_gate(c) & not_gate(d));
                lo=lx1;
                lg=lx2;
                break;

            case LUT_AXBM_EN_2:
                //a=di+2 b=di+1 c=d d=di-1
                lx1=(not_gate(a) & not_gate(b) & not_gate(c)) | (a & b & d) | (not_gate(a) & not_gate(b) & c & not_gate(d))|(a & b & c & not_gate(d));
                lx2=(not_gate(b) & c & d) | (not_gate(a) & b & not_gate(c) & not_gate(d)) | (a & b & not_gate(c) & not_gate(d))| (a & not_gate(b) & not_gate(c) & d) | (a & not_gate(b) & c & not_gate(d));
                lo=lx1;
                lg=lx2;
                break;
            case LUT_BW_A:
                //f=a_i  e=b_i d=pp
                sig1=f & e;
                sig2=d ^ (sig1);
                lo=sig2;
                lg=d;
                break;
            case LUT_BW_B:
                //f=a_i  e=b_i d=pp c=0/1
                sig1=not_gate(f & e);
                sig2=d ^ c;
                lo= sig1 ^ sig2;
                lg=sig2;
                break;
            case LUT_BW_C:
                lo=1;
                lg=1;
                break;
            case LUT_BWO_J:
                //f=a_i  e=b_i d=pp
                sig1=f & e; //Upper part
                sig2=c & d; //Lower part
                lo=b ^ sig1 ^ sig2;
                //lg=b ^ sig1;
                lg=b;
                break;
            case LUT_BWO_K:
                //f=a_i  e=b_i d=pp c=0/1
                sig1=not_gate(f & e); //Upper part
                sig2=c & d;  //Lower part
                lo= b ^ sig1 ^ sig2;
                //lg=b ^ sig1;
                lg=b;
                break;
            case LUT_BWO_L:
                //f=a_i  e=b_i d=pp c=0/1
                sig1 = f & e; //Upper part
                sig2 = not_gate(c & d);  //Lower part
                lo = b ^ sig1 ^ sig2;
                //lg = b ^ sig1;
                lg=b;
                break;
             case LUT_BWO_M:
                //f=a_i  e=b_i d=pp c=0/1
                sig1=not_gate(f & e);
                sig2=not_gate(c & d);
                lo= b ^ sig1 ^ sig2;
                //lg=b ^ sig1;
                lg=b;
                break;
            case LUT_BWO_N:
                lo=1;
                lg=1;
                break;
            case LUT_TEST:
                lo=not_gate(a) & not_gate(b) & not_gate(c) & not_gate(d) & not_gate(e) & f;
                break;
        }
        *counter=*counter+1;

        result[0]=lo;
        result[1]=lg;
}


void SuperLUT_6_2(uchar a, uchar b, uchar c, uchar d, uchar e, uchar f, uchar lut_type, int* counter, uchar cin,uchar* cout,uchar* result) //2 equiprobale LUTs
{

  uchar luti_result[2]={ };
  uchar luts_result[2]={ };
  LUT_6_2(a,b,c,d,e,f,lut_type,counter,luti_result);
  uchar luts_cin=cin;
  uchar res=0;
  for(uchar cc=0;cc<2;cc++)
  {
    CARRY1(luti_result[0],luti_result[1],luts_cin,&luts_cin,&res);
    luts_result[cc]=res;
  }
    result[0]=luts_result[0];
    result[1]=luts_result[1];
    *cout=luts_cin;
}


void lut_testing()
{
    int lut_counter=0;
    uchar o[2]={ };
    for(uchar a=0;a<=1;a++)
        for(uchar b=0;b<=1;b++)
            for(uchar c=0;c<=1;c++)
                for(uchar d=0;d<=1;d++)
                {
                    LUT_6_2(a,b,c,d,0,0,LUT_AXBM_EN_1,&lut_counter,o);
                    printf("a=%d b=%d c=%d d=%d x1=%d x2=%d \n",a,b,c,d,o[0],o[1]);
                }
}


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

long efficient_softcore(long multiplicand, long multiplier)
{
    long multiplier_appended=multiplier<<1;
    long multiplicand_appended=multiplicand<<2;
    uchar bits[3]={ };
    long pp_products[BIT_WIDTH*2+4]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[BIT_WIDTH+4]={ }, t_row1[BIT_WIDTH+4]={ };//initialize all zero
    t_row[BIT_WIDTH+3]=1;
    t_row[BIT_WIDTH+2]=1;
    uchar lut_result[2]={0,0};
    uchar cin=1,res=0,pp_in=0;
    long prod=0;
    for(int i=0;i<PP_NUM_UNSIGNED;i++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended>>(2*i) & 1;
        bits[BIT_CSB]=multiplier_appended>>(2*i+1) & 1;
        bits[BIT_MSB]=multiplier_appended>>(2*i+2) & 1;
        cin=1;
        for(int pp=0;pp<BIT_WIDTH+4;pp++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp+1) & 1;
            pp_in=t_row[pp];
            if(pp==(BIT_WIDTH+3))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_ES_C,&lut_counter,lut_result);
            }
            else if(pp==(BIT_WIDTH+2))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_ES_B,&lut_counter,lut_result);
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_ES_A,&lut_counter,lut_result);
            }
            //goes in to carry
            CARRY1(lut_result[0],lut_result[1],cin,&cin,&res);
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp<3 && pp>0)
            {
                pp_products[i*2+pp-1]=res;
            }else if(pp>=3)
            {
                t_row1[pp-2]=res; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }
        t_row1[BIT_WIDTH+2]=cin;
        for(int t=0; t<BIT_WIDTH+4;t++){t_row[t]=t_row1[t];}
    }

    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_UNSIGNED*2)
            {
                    prod|=(long)pp_products[i]<<i;
            }
            else{
                    prod|=(long)t_row[i-(PP_NUM_UNSIGNED*2)+1] << i;
            }
    }
    return(prod);
}


long efficient_softcore_hdl(long multiplicand, long multiplier)
{
    long multiplier_appended=multiplier<<1;
    long multiplicand_appended=multiplicand<<2;
    uchar bits[3]={ };
    long pp_products[PP_NUM_UNSIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_UNSIGNED+1][BIT_WIDTH+4]={ };//initialize all zero
    t_row[0][BIT_WIDTH+3]=1;
    t_row[0][BIT_WIDTH+2]=1;
    uchar lut_result[PP_NUM_UNSIGNED][BIT_WIDTH+4][2]={0,0};
    uchar carry[PP_NUM_UNSIGNED][BIT_WIDTH+5]={ };
    uchar res[PP_NUM_UNSIGNED][BIT_WIDTH+4]={ };
    long prod=0;
    for(int pp_row=0;pp_row<PP_NUM_UNSIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        carry[pp_row][0]=1;
        for(int pp_col=0;pp_col<BIT_WIDTH+4;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            if(pp_col==(BIT_WIDTH+3))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_ES_C,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH+2))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_ES_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_ES_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col],&carry[pp_row][pp_col+1],&res[pp_row][pp_col]);
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<3 && pp_col>0)
            {
                pp_products[pp_row*2+pp_col-1]=res[pp_row][pp_col];
            }else if(pp_col>=3)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+2]=carry[pp_row][BIT_WIDTH+4];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_UNSIGNED*2)
            {
                    prod|=(long)pp_products[i]<<i;
            }
            else{
                    prod|=(long)t_row[PP_NUM_UNSIGNED][i-(PP_NUM_UNSIGNED*2)+1] << i;
            }
    }
    return(prod);
}



long long baugh_wooley_hdl(long long multiplicand, long long multiplier)
{
    uchar bits[3]={ };
    long pp_products[BIT_WIDTH]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[BIT_WIDTH+1][BIT_WIDTH+1]={ };//initialize all zero

    uchar lut_result[BIT_WIDTH][BIT_WIDTH+1][2]={0,0};
    //uchar carry[PP_NUM_UNSIGNED][BIT_WIDTH+5]={ };
    uchar carry=0;
    uchar res[BIT_WIDTH][BIT_WIDTH+1]={ };
    long prod=0;
    for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        carry=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
        {
            multiple_tuple[0]=(multiplicand >> pp_col) & 1;
            multiple_tuple[1]=(multiplier >> pp_row) & 1;
            //######1
            if(pp_row == BIT_WIDTH-1)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    //A
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else if(pp_col==BIT_WIDTH)
                {
                    //C
                    LUT_6_2(0,0,0,0,0,0,LUT_BW_C,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                    //B c->0
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
            }
           //########2
            if(pp_row==1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    //B c->1
                    LUT_6_2(0,0,1,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                    //A
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                }
            }
            //#########3
            else if(pp_row<BIT_WIDTH-1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                        //B c->0
                        LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                        //A
                        LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                }
            }
            //goes in to carry
            if((pp_row < (BIT_WIDTH-1) && pp_col < (BIT_WIDTH)) || (pp_row == (BIT_WIDTH-1)))
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry,&carry,&res[pp_row][pp_col]);
            }
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col==0)
            {
                pp_products[pp_row]=res[pp_row][pp_col];
            }else
            {
                t_row[pp_row+1][pp_col-1]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop

        if(pp_row==BIT_WIDTH-1)
        {
            t_row[pp_row+1][BIT_WIDTH-1]=res[pp_row][BIT_WIDTH];
        }
        else{
            t_row[pp_row+1][BIT_WIDTH-1]=carry;
        }
    }//pp_row loop
    //Last row:


    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<BIT_WIDTH)
            {
                    prod|=(long)pp_products[i]<<i;
            }
            else{
                    prod|=(long)t_row[BIT_WIDTH][i-BIT_WIDTH] << i;
            }
    }
#ifdef CARRY_PROBE
        for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
        {
            for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
            {
                pp_probe[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
        }
        //##probe pp_product
#endif // CARRY_PROBE

    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}

long long lc_baugh_wooley_opt(long long multiplicand, long long multiplier, int version)
{
    uchar bits[3]={ };
    long pp_products[BIT_WIDTH]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[BIT_WIDTH+1][BIT_WIDTH+1]={ };//initialize all zero

    uchar lut_result[BIT_WIDTH][BIT_WIDTH+1][2]={0,0};
    //uchar carry[PP_NUM_UNSIGNED][BIT_WIDTH+5]={ };
    uchar carry=0;
    uchar res[BIT_WIDTH][BIT_WIDTH+1]={ };
    long prod=0;
    int cnt=0;
    uchar lut_res[2]={ };
    int ver=version+1;
    for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        carry=0;
        cnt=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
        {
            multiple_tuple[0]=(multiplicand >> pp_col) & 1;
            multiple_tuple[1]=(multiplier >> pp_row) & 1;

            if(pp_row == BIT_WIDTH-1)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else if(pp_col==BIT_WIDTH)
                {
                    //C
                    LUT_6_2(0,0,0,0,0,0,LUT_BW_C,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                    //B c->0
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }


            }
            if(pp_row==1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    //B c->1
                    LUT_6_2(0,0,1,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                        if(pp_col<=(BIT_WIDTH-ver-pp_row))
                        {
                         if(cnt==0)
                         {
                            multiple_tuple[0]=(multiplicand >> BIT_WIDTH-ver-pp_row) & 1;
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_res);
                         }
                         cnt++;
                         lut_result[pp_row][pp_col][0]=lut_res[0];
                         lut_result[pp_row][pp_col][1]=lut_res[1];
                        }
                        else{
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                        }


                }


            }
            else if(pp_row<BIT_WIDTH-1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                        //B c->0
                        LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                        //A
                        if(pp_col<=(BIT_WIDTH-ver-pp_row))
                        {
                         if(cnt==0)
                         {
                            multiple_tuple[0]=(multiplicand >> BIT_WIDTH-ver-pp_row) & 1;
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_res);
                         }
                         cnt++;
                         lut_result[pp_row][pp_col][0]=lut_res[0];
                         lut_result[pp_row][pp_col][1]=lut_res[1];
                        }
                        else{
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                        }
                }
            }

            //goes in to carry
            if((pp_row < (BIT_WIDTH-1) && pp_col < (BIT_WIDTH)) || (pp_row == (BIT_WIDTH-1)))
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry,&carry,&res[pp_row][pp_col]);
            }
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col==0)
            {
                pp_products[pp_row]=res[pp_row][pp_col];
            }else
            {
                t_row[pp_row+1][pp_col-1]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop

        if(pp_row==BIT_WIDTH-1)
        {
            t_row[pp_row+1][BIT_WIDTH-1]=res[pp_row][BIT_WIDTH];
        }
        else{
            t_row[pp_row+1][BIT_WIDTH-1]=carry;
        }
    }//pp_row loop
    //Last row:


    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<BIT_WIDTH)
            {
                    prod|=(long)pp_products[i]<<i;
            }
            else{
                    prod|=(long)t_row[BIT_WIDTH][i-BIT_WIDTH] << i;
            }
    }
#ifdef CARRY_PROBE
        for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
        {
            for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
            {
                pp_probe[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
        }
        //##probe pp_product
#endif // CARRY_PROBE

    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}

long long lc_baugh_wooley_opt_ac(long long multiplicand, long long multiplier, int version)
{
    uchar bits[3]={ };
    long pp_products[BIT_WIDTH]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[BIT_WIDTH+1][BIT_WIDTH+1]={ };//initialize all zero

    uchar lut_result[BIT_WIDTH][BIT_WIDTH+1][2]={0,0};
    //uchar carry[PP_NUM_UNSIGNED][BIT_WIDTH+5]={ };
    uchar carry=0;
    uchar res[BIT_WIDTH][BIT_WIDTH+1]={ };
    long prod=0;
    int cnt=0;
    uchar lut_res[2]={ };
    int ver=version+1;
    uchar carry_temp=0;
    uchar res_temp=0;
    for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        carry=0;
        cnt=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
        {
            multiple_tuple[0]=(multiplicand >> pp_col) & 1;
            multiple_tuple[1]=(multiplier >> pp_row) & 1;

            if(pp_row == BIT_WIDTH-1)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else if(pp_col==BIT_WIDTH)
                {
                    //C
                    LUT_6_2(0,0,0,0,0,0,LUT_BW_C,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                    //B c->0
                    LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }


            }
            if(pp_row==1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                    //B c->1
                    LUT_6_2(0,0,1,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                        if(pp_col<=(BIT_WIDTH-ver-pp_row))
                        {
                         if(cnt==0)
                         {
                            multiple_tuple[0]=(multiplicand >> BIT_WIDTH-ver-pp_row) & 1;
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_res);
                         }
                         //cnt++;
                         lut_result[pp_row][pp_col][0]=lut_res[0];
                         lut_result[pp_row][pp_col][1]=lut_res[1];
                        }
                        else{
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                        }


                }


            }
            else if(pp_row<BIT_WIDTH-1 && pp_col < BIT_WIDTH)
            {
                if(pp_col==BIT_WIDTH-1)
                {
                        //B c->0
                        LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_B,&lut_counter,lut_result[pp_row][pp_col]);
                }
                else
                {
                        //A
                        if(pp_col<=(BIT_WIDTH-ver-pp_row))
                        {
                         if(cnt==0)
                         {
                            multiple_tuple[0]=(multiplicand >> BIT_WIDTH-ver-pp_row) & 1;
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_res);
                         }
                         //cnt++;
                         lut_result[pp_row][pp_col][0]=lut_res[0];
                         lut_result[pp_row][pp_col][1]=lut_res[1];
                        }
                        else{
                            LUT_6_2(0,0,0,t_row[pp_row][pp_col],multiple_tuple[1],multiple_tuple[0],LUT_BW_A,&lut_counter,lut_result[pp_row][pp_col]);
                        }
                }
            }

            //goes in to carry
            if((pp_row < (BIT_WIDTH-1) && pp_col < (BIT_WIDTH)) || (pp_row == (BIT_WIDTH-1)))
            {
              if(pp_col<=(BIT_WIDTH-ver-pp_row))
              {
                if(cnt==0)
                {
                  CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry,&carry_temp,&res_temp);
                }
                cnt++;
                carry=carry_temp;
                res[pp_row][pp_col]=res_temp;
             }
             else{
                CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry,&carry,&res[pp_row][pp_col]);
             }

            }
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col==0)
            {
                pp_products[pp_row]=res[pp_row][pp_col];
            }else
            {
                t_row[pp_row+1][pp_col-1]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop

        if(pp_row==BIT_WIDTH-1)
        {
            t_row[pp_row+1][BIT_WIDTH-1]=res[pp_row][BIT_WIDTH];
        }
        else{
            t_row[pp_row+1][BIT_WIDTH-1]=carry;
        }
    }//pp_row loop
    //Last row:


    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<BIT_WIDTH)
            {
                    prod|=(long)pp_products[i]<<i;
            }
            else{
                    prod|=(long)t_row[BIT_WIDTH][i-BIT_WIDTH] << i;
            }
    }
#ifdef CARRY_PROBE
        for(int pp_row=0;pp_row<BIT_WIDTH;pp_row++)
        {
            for(int pp_col=0;pp_col<BIT_WIDTH+1;pp_col++)
            {
                pp_probe[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
        }
        //##probe pp_product
#endif // CARRY_PROBE

    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}





//signed

long long area_optimized_hdl(long long multiplicand, long long multiplier)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<2;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+4]={ };//initialize all zero
    for(int i=0;i<PP_NUM_SIGNED+1;i++){t_row[i][BIT_WIDTH+3]=1;}
    t_row[0][BIT_WIDTH+2]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+4][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+5]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+4]={ };
    long long prod=0;
    uchar pp_in=0;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        carry[pp_row][0]=1; //this has to be intialized for all the rows in vhdl
        for(int pp_col=0;pp_col<BIT_WIDTH+4;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==(BIT_WIDTH+3))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_Z,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH+2))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col],&carry[pp_row][pp_col+1],&res[pp_row][pp_col]);
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<3 && pp_col>0)
            {
                pp_products[pp_row*2+pp_col-1]=res[pp_row][pp_col];
            }else if(pp_col>=3)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+2]=carry[pp_row][BIT_WIDTH+4];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)+1] << i;
            }
    }
    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}



// The B LUT for last PP row is not required hence it can be removed.
long long area_optimized_hdl_opt(long long multiplicand, long long multiplier)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    t_row[0][BIT_WIDTH+1]=1;
    t_row[0][BIT_WIDTH]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+2][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    long long prod=0;
    uchar pp_in=0;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        for(int pp_col=0;pp_col<BIT_WIDTH+2;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==0)
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
             carry[pp_row][pp_col]=lut_result[pp_row][pp_col][1];
             res[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
            else if(pp_col==(BIT_WIDTH+1))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            if(pp_col>0)
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry[pp_row][pp_col],&res[pp_row][pp_col]);
            }

            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res[pp_row][pp_col];
            }else if(pp_col>=2)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+1]=carry[pp_row][BIT_WIDTH+1];
        t_row[pp_row+1][BIT_WIDTH]=carry[pp_row][BIT_WIDTH+1];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }

#ifdef CARRY_PROBE

    //Carry probe

   for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_opt[i1][j1] = carry[i1][j1] ;
        }

        }

        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_opt[i1][j1] = lut_result[i1][j1][0];
        }

        }
#endif
    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}


long long area_optimized_hdl_approx(long long multiplicand, long long multiplier)
{

    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    t_row[0][BIT_WIDTH+1]=1;
    t_row[0][BIT_WIDTH]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+2][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    long long prod=0;
    uchar pp_in=0;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        for(int pp_col=0;pp_col<BIT_WIDTH+2;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==0)//
            {
             if(pp_row==PP_NUM_SIGNED-1)
             {
                LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
                carry[pp_row][pp_col]=lut_result[pp_row][pp_col][1];
                res[pp_row][pp_col]=0;
             }
             else
             {
                carry[pp_row][pp_col]=1;
                res[pp_row][pp_col]=0;
             }
            }
            else if(pp_col==1)//
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_Am,&lut_counter,lut_result[pp_row][pp_col]);
             carry[pp_row][pp_col]=carry[pp_row][pp_col-1];
             res[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
            else if(pp_col==(BIT_WIDTH+1))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            if(pp_col>1)
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry[pp_row][pp_col],&res[pp_row][pp_col]);
            }

            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res[pp_row][pp_col];
            }else if(pp_col>=2)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+1]=carry[pp_row][BIT_WIDTH+1];
        t_row[pp_row+1][BIT_WIDTH]=carry[pp_row][BIT_WIDTH+1];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }
    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}

long long lc_booth_opt(long long multiplicand, long long multiplier,int version)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    t_row[0][BIT_WIDTH+1]=1;
    t_row[0][BIT_WIDTH]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+2][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    long long prod=0;
    uchar pp_in=0;
    uchar lut_temp[2]={ };
    int cnt=0;
    int ver=2-version;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        cnt=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+2;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==0)
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
             carry[pp_row][pp_col]=lut_result[pp_row][pp_col][1];
             res[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
            else if(pp_col==(BIT_WIDTH+1))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col>0 && pp_col<=((2*(PP_NUM_SIGNED-pp_row-1)+ver)))
            {
                if(cnt==0)
                {
                  multiple_tuple[0]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+ver) & 1);
                  multiple_tuple[1]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+ver+1) & 1);
                  LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_temp);
                    //LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,1,0,LUT_AO_A,&lut_counter,lut_temp);
                }
                cnt++;
                lut_result[pp_row][pp_col][0]=lut_temp[0];
                lut_result[pp_row][pp_col][1]=lut_temp[1];
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            if(pp_col>0)
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry[pp_row][pp_col],&res[pp_row][pp_col]);
            }

            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res[pp_row][pp_col];
            }else if(pp_col>=2)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+1]=carry[pp_row][BIT_WIDTH+1];
        t_row[pp_row+1][BIT_WIDTH]=carry[pp_row][BIT_WIDTH+1];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }

#ifdef CARRY_PROBE

    //Carry probe

   for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_opt[i1][j1] = carry[i1][j1] ;
        }

        }

        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_opt[i1][j1] = lut_result[i1][j1][0] ;
        }

        }
#endif
    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}



long long lc_booth_opt_ac(long long multiplicand, long long multiplier,int version)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    t_row[0][BIT_WIDTH+1]=1;
    t_row[0][BIT_WIDTH]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+2][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    long long prod=0;
    uchar pp_in=0;
    uchar lut_temp[2]={ };
    int cnt=0;
    int ver=2-version;
    uchar carry_temp=0;
    uchar res_temp=0;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        cnt=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+2;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==0)
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
             carry[pp_row][pp_col]=lut_result[pp_row][pp_col][1];
             res[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
            else if(pp_col==(BIT_WIDTH+1))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col>0 && pp_col<=((2*(PP_NUM_SIGNED-pp_row-1)+ver)))
            {
                if(cnt==0)
                {
                  multiple_tuple[0]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+ver) & 1);
                  multiple_tuple[1]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+ver+1) & 1);
                  LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_temp);
                    //LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,1,0,LUT_AO_A,&lut_counter,lut_temp);
                }
                lut_result[pp_row][pp_col][0]=lut_temp[0];
                lut_result[pp_row][pp_col][1]=lut_temp[1];
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            if(pp_col>0)
            {
            //################ In test ##################
             if(pp_col<=((2*(PP_NUM_SIGNED-pp_row-1)+ver)))
             {
              if(cnt==0)
              {
                CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry_temp,&res_temp);
              }
                cnt++;
                carry[pp_row][pp_col]=carry_temp;
                res[pp_row][pp_col]=res_temp;
             }
             else{
                CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry[pp_row][pp_col],&res[pp_row][pp_col]);
             }
             //##########################################

            }

            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res[pp_row][pp_col];
            }else if(pp_col>=2)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+1]=carry[pp_row][BIT_WIDTH+1];
        t_row[pp_row+1][BIT_WIDTH]=carry[pp_row][BIT_WIDTH+1];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }

#ifdef CARRY_PROBE

    //Carry probe

   for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_opt[i1][j1] = carry[i1][j1] ;
        }

        }

        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_opt[i1][j1] = lut_result[i1][j1][0] ;
        }

        }
#endif
    prod=unsigned_to_signed(prod,2*BIT_WIDTH);
    return(prod);
}






/*
long long stair_booth_opt_v2_2(long long multiplicand, long long multiplier)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    t_row[0][BIT_WIDTH+1]=1;
    t_row[0][BIT_WIDTH]=1;
    uchar lut_result[PP_NUM_SIGNED][BIT_WIDTH+2][2]={0,0};
    uchar carry[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    uchar res[PP_NUM_SIGNED][BIT_WIDTH+2]={ };
    long long prod=0;
    uchar pp_in=0;
    uchar lut_temp[2]={ };
    int cnt=0;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        cnt=0;
        for(int pp_col=0;pp_col<BIT_WIDTH+2;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> pp_col) & 1;
            multiple_tuple[1]=(multiplicand_appended >> pp_col+1) & 1;
            pp_in=t_row[pp_row][pp_col];
            if(pp_col==0)
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
             carry[pp_row][pp_col]=lut_result[pp_row][pp_col][1];
             res[pp_row][pp_col]=lut_result[pp_row][pp_col][0];
            }
            else if(pp_col==(BIT_WIDTH+1))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(BIT_WIDTH))
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col>0 && pp_col<=((2*(PP_NUM_SIGNED-pp_row-1)+1)))
            {
                if(cnt==0)
                {
                  multiple_tuple[0]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+1) & 1);
                  multiple_tuple[1]=(multiplicand_appended >> (2*(PP_NUM_SIGNED-pp_row-1)+2) & 1);
                  LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_temp);
                    //LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,1,0,LUT_AO_A,&lut_counter,lut_temp);
                }
                cnt++;
                lut_result[pp_row][pp_col][0]=lut_temp[0];
                lut_result[pp_row][pp_col][1]=lut_temp[1];
            }
            else
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }
            //goes in to carry
            if(pp_col>0)
            {
             CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],carry[pp_row][pp_col-1],&carry[pp_row][pp_col],&res[pp_row][pp_col]);
            }

            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res[pp_row][pp_col];
            }else if(pp_col>=2)
            {
                t_row[pp_row+1][pp_col-2]=res[pp_row][pp_col]; //the t_row[0] is supposed to be zero (LSB LUT of PP)
            }
        }//pp_col loop
        t_row[pp_row+1][BIT_WIDTH+1]=carry[pp_row][BIT_WIDTH+1];
        t_row[pp_row+1][BIT_WIDTH]=carry[pp_row][BIT_WIDTH+1];
    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }

#ifdef CARRY_PROBE

    //Carry probe

   for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_opt[i1][j1] = carry[i1][j1] ;
        }

        }

        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_opt[i1][j1] = lut_result[i1][j1][0] ;
        }

        }
#endif

    return(prod);
}
*/
/*
long long star_booth_opt_v3(long long multiplicand, long long multiplier)
{
    long long multiplier_appended=multiplier<<1;
    long long multiplicand_appended=multiplicand<<1;
    uchar bits[3]={ };
    long long pp_products[PP_NUM_SIGNED*2]={ };
    int lut_counter=0;
    uchar multiple_tuple[2];
    //condition t_row
    //uchar T_NUM=
    //WHen implementing VHDL this can be precisely calculated.
    uchar t_row[PP_NUM_SIGNED+1][BIT_WIDTH+2]={ };//initialize all zero
    uchar lut_result[PP_NUM_SIGNED+1][BIT_WIDTH+2][2]={ };
    //uchar lut
    long long prod=0;
    uchar pp_in=0;
    uchar lut_temp[2]={ };
    uchar pp_row_cols=5; //
    uchar pp_slab=0;
    uchar cin=0;
    uchar sulut1_result[3]={ };
    uchar sulut2_result[2]={ };
    uchar res=0;
    t_row[0][4]=1;
    t_row[0][3]=1;
    for(int pp_row=0;pp_row<PP_NUM_SIGNED;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(2*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(2*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(2*pp_row+2) & 1;
        pp_slab=0;
        cin=0;
        for(int pp_col=0;pp_col<pp_row_cols;pp_col++)
        {
            multiple_tuple[0]=(multiplicand_appended >> (pp_col+pp_slab)) & 1;  //SA lut multiplie tuple[0] : pp_col value=1 + (slab-1)
            multiple_tuple[1]=(multiplicand_appended >> (pp_col+pp_slab+1)) & 1;    //SA lut multiplie tuple[0] : pp_col value=1 + (slab)
            pp_in=t_row[pp_row][pp_col]; //Have to adjust this.
            if(pp_col==0)
            {
             LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A1,&lut_counter,lut_result[pp_row][pp_col]);
             cin=lut_result[pp_row][pp_col][1];
             res=lut_result[pp_row][pp_col][0];
             if(pp_row < (PP_NUM_SIGNED-1))
             {
                //pp_slab=(2*(PP_NUM_SIGNED-pp_row-1))-1;

             }
            }
            else if(pp_row < (PP_NUM_SIGNED-2) && pp_col==1)
            {
                //SuperLUT1_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,cin,&cin,sulut1_result,pp_slab);
                SuperLUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,cin,&cin,sulut2_result);
                res=sulut2_result[0];
                t_row[pp_row+1][0]=sulut2_result[1];
                t_row[pp_row+1][1]=sulut2_result[1];
            }
            else if(pp_row == (PP_NUM_SIGNED-2) && pp_col==1)
            {
                SuperLUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,cin,&cin,sulut2_result);
                res=sulut2_result[0];
                t_row[pp_row+1][0]=sulut2_result[1];
            }
            else if(pp_col==(pp_row_cols-1))
            {
                LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_B,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else if(pp_col==(pp_row_cols-2))
            {
                LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A2,&lut_counter,lut_result[pp_row][pp_col]);
            }
            else
            {
                LUT_6_2(bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],pp_in,multiple_tuple[1],multiple_tuple[0],LUT_AO_A,&lut_counter,lut_result[pp_row][pp_col]);
            }

            //goes in to carry
            if(pp_row==PP_NUM_SIGNED-1)
            {
                if(pp_col>0) //0->A1 1->SA1/2 No carry block required.
                {
                    CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],cin,&cin,&res);
                }

            }
            else{
                if(pp_col>1) //0->A1 1->SA1/2 No carry block required.
                {
                    CARRY1(lut_result[pp_row][pp_col][0],lut_result[pp_row][pp_col][1],cin,&cin,&res);
                }
            }
            //accumulate results x,p1,p2,t0,t1,........... t9,1
            if(pp_col<2)
            {
                pp_products[pp_row*2+pp_col]=res;
            }else if(pp_col>=2)
            {
                if(pp_row==PP_NUM_SIGNED-2)
                    t_row[pp_row+1][pp_col-1]=res; //the t_row[0] is supposed to be zero (LSB LUT of PP)
                else if(pp_row==PP_NUM_SIGNED-1)
                    t_row[pp_row+1][pp_col-2]=res;
                else
                    t_row[pp_row+1][pp_col]=res;
            }

        }//pp_col loop

        if(pp_row==(PP_NUM_SIGNED-2))
        {
            t_row[pp_row+1][pp_row_cols-1]=cin;
            t_row[pp_row+1][pp_row_cols]=cin;
        }
        else
        {
            t_row[pp_row+1][pp_row_cols]=cin;
            t_row[pp_row+1][pp_row_cols+1]=cin;
        }

        if(pp_row<(PP_NUM_SIGNED-2))
        {
            pp_row_cols+=2;
        }
        else{
            pp_row_cols+=1;
        }

    }//pp_row loop
    for(int i=0;i<BIT_WIDTH*2;i++)
    {
            if(i<PP_NUM_SIGNED*2)
            {
                    prod|=(long long)pp_products[i]<<i;
            }
            else{
                    prod|=(long long)t_row[PP_NUM_SIGNED][i-(PP_NUM_SIGNED*2)] << i;
            }
    }

    /*
#ifdef CARRY_PROBE

    //Carry probe

   for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_opt[i1][j1] = carry[i1][j1] ;
        }

        }

        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_opt[i1][j1] = lut_result[i1][j1][0] ;
        }

        }
#endif

    return(prod);
}
*/



void probe_carry(uchar *carry_values)
{
        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            carry_values[i1*(BIT_WIDTH+2)+j1] = carry_opt[i1][j1];
        }
        }
}

void probe_lut_result(uchar *lut_values)
{
        for (int i1 = 0; i1 < PP_NUM_SIGNED; i1++){
        for (int j1 = 0; j1 < (BIT_WIDTH+2); j1++){
            lut_values[i1*(BIT_WIDTH+2)+j1] = lut_opt[i1][j1];
        }
        }
}


void probe_values(uchar *lut_values)
{
        for (int i1 = 0; i1 < BIT_WIDTH; i1++){
        for (int j1 = 0; j1 < 2*BIT_WIDTH; j1++){
            lut_values[i1*2*BIT_WIDTH+j1] = pp_probe[i1][j1];
        }
        }
}




////////////////////////IEEE Access: FPGA-Based Multi-Level Approximate Multipliers for High-Performance Error-Resilient Applications

void half_adder(uchar a, uchar b, uchar* sum, uchar* cout)
{
	uchar result;
	uchar ccout;
	result = a ^ b;
	ccout = a & b;
    *cout=ccout;
	*sum=result;
}

void full_adder(uchar a, uchar b, uchar cin, uchar* sum, uchar* cout)
{
	uchar result;
	uchar ccout;
	result = a ^ b ^ cin;
	ccout = (a & b) ^( cin & (a ^ b));
    *cout=ccout;
	*sum=result;
}


void exact_compressor(uchar a, uchar b, uchar c, uchar d, uchar cin, uchar *sum, uchar *cout, uchar *carry)
{
    uchar co=0;
    uchar carr=0;
    uchar s=0;
    uchar o1=0;
    full_adder(a, b, c, &co, &o1);
    full_adder(o1, d, cin, &carr, &s);
    *sum=s;
    *cout=co;
    *carry=carr;
}


// GOLDEN FOR SIGNED
void approximate_compressor_3_2(uchar a, uchar b, uchar c, uchar* res1, uchar* res2)
{
    //swapper a<->c
    a=a+c;c=a-c;a=a-c;
    uchar result[2]={ };
    result[0]= a | b;
    result[1]=c;

    *res1=result[0];
    *res2=result[1];
}



void approximate_compressor_4_2_c_1(uchar a, uchar b, uchar c, uchar d, uchar* res1, uchar* res2)
{
    //swapper a<->d ,b<->c
    a=a+d;d=a-d;a=a-d;
    b=b+c;c=b-c;b=b-c;
    //verified
    uchar result[2]={ };
    result[0]=((d | c) & (b | a)) | (d & c) | (b & a); //C
    result[1]=(a ^ b ^ c ^ d) | (a & b & c & d); //S

    *res1=result[0];
    *res2=result[1];
}

void approximate_compressor_4_2_c_2(uchar a, uchar b, uchar c, uchar d, uchar* res1, uchar* res2)
{
    //swapper a<->d ,b<->c
    //a=a+d;d=a-d;a=a-d;
    //b=b+c;c=b-c;b=b-c;


    uchar result[2]={ };
    result[0]=((d | c) & (b | a)) | (d & c);
    result[1]=(b | a) | (d ^ c);

    *res1=result[0];
    *res2=result[1];
}



void approximate_compressor_4_2_c_3(uchar a, uchar b, uchar c, uchar d, uchar* res1, uchar* res2)
{
    //swapper a<->d ,b<->c
    a=a+d;d=a-d;a=a-d;
    b=b+c;c=b-c;b=b-c;


    uchar result[2]={ };
    result[0]=b | a;
    result[1]=d | c;


    *res1=result[0];
    *res2=result[1];
}


void approximate_compressor_4_2_c_4(uchar a, uchar b, uchar c, uchar d, uchar* res1, uchar* res2)
{
    //swapper a<->d ,b<->c
    a=a+d;d=a-d;a=a-d;
    b=b+c;c=b-c;b=b-c;

    uchar result[2]={ };
    result[0]=b | a;
    result[1]=d | c | (b & a);

    *res1=result[0];
    *res2=result[1];
}



long long multilevel_signed_opt(long long multiplicand_sign, long long multiplier_sign, char c2, char k_paper)
{
    long long factor=(long long)(pow(2,BIT_WIDTH_ML)-1);
    long long multiplicand=multiplicand_sign & factor;
    long long multiplier=multiplier_sign & factor;

    char k=k_paper+2;
    uchar row_stop=0;
    long long pps_acc[BIT_WIDTH_ML]={ };
    long long mask=0;
    long long pps[BIT_WIDTH_ML]={ };
    long long pps_r=0;
    long long pp_product=0;
    long long pps_st2[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    long long pps_st3[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    long long pps_st1[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    uchar carry[BIT_WIDTH_ML/4][2*BIT_WIDTH_ML]={ };
    uchar pp_n[2*BIT_WIDTH_ML]={ };
    uchar pp_n_l[2*BIT_WIDTH_ML]={ };
    int cnt=0;
    uchar pp_out=0;
    int pp_row1=0;
    uchar cin=0;
    long long multiplicand_in=multiplicand;
    long long multiplier_in=multiplier;
    //accurate part: NOTE it also requires the carry from the previous part
    char r_shift=0;
    for(int pp_row=0; pp_row<BIT_WIDTH_ML ;pp_row++)
    {
        pps_r=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        //complement the leading bits using XOR gate

        if(pp_row==BIT_WIDTH_ML-1)
        {
            mask=(1 << (BIT_WIDTH_ML-1));
            mask=(mask-1) |  (1 << (BIT_WIDTH_ML));
            pps_r=pps_r ^ mask;
        }
        else if(pp_row==0)
        {
            pps_r=pps_r ^ (3 << (BIT_WIDTH_ML-1));
        }
        else
        {
            pps_r=pps_r ^ (1 << (BIT_WIDTH_ML-1));
        }
//        pps_acc[pp_row]=pps_r<<pp_row; //ONLY for test
        r_shift=(2*BIT_WIDTH_ML)-k-pp_row;
        if(r_shift<0)
        {
            pps_acc[pp_row]=pps_r<<((2*BIT_WIDTH_ML)-k-r_shift);
        }
        else{
                pps_acc[pp_row]=(pps_r>>r_shift)<<((2*BIT_WIDTH_ML)-k);
        }

        //print_binary2(pps_acc[pp_row],16);
        //printf("\n");

        pp_product+=pps_acc[pp_row];
    }

    //Have to shift carry fromt he approximate part by << (2*BIT_WIDTH_ML)-k;
    pp_product&=(1 << 2*BIT_WIDTH_ML)-1;
   // print_binary2(pp_product,16);
   // printf("\n");
   // return(pp_product);

    //approximate part
    for(int i=1;i<=BIT_WIDTH_ML;i++)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int i=BIT_WIDTH_ML-1;i>0;i--)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int pp_row=0; pp_row<BIT_WIDTH_ML ;pp_row++)
    {
        pps[pp_row]=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        if(pp_row==BIT_WIDTH_ML-1)
        {
            mask=(1 << (BIT_WIDTH_ML-1));
            mask=(mask-1) |  (1 << (BIT_WIDTH_ML));
            pps[pp_row]=pps[pp_row] ^ mask;
        }
        else if(pp_row==0)
        {
            pps[pp_row]=pps[pp_row] ^ (3 << (BIT_WIDTH_ML-1));
        }
        else
        {
            pps[pp_row]=pps[pp_row] ^ (1 << (BIT_WIDTH_ML-1));
        }
        pp_row1=pp_row;
        for(int pp_col=pp_row; pp_col < BIT_WIDTH_ML+pp_row;pp_col++)
        {
            if(pp_col>=BIT_WIDTH_ML)
            {
                pps_st1[--pp_row1][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
                //pp_row1--;
            }else
            {
                pps_st1[pp_row][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
            }

        }

    }

    //insert 1 at he BIT_WIDTH_ML+1 position ALSO increase the counter at this position.
    pp_n[BIT_WIDTH_ML]+=1;
    for(int pp_row=BIT_WIDTH_ML-1;pp_row>0;pp_row--)
    {
        pps_st1[pp_row][BIT_WIDTH_ML]=pps_st1[pp_row-1][BIT_WIDTH_ML];
    }
    pps_st1[0][BIT_WIDTH_ML]=1;
    pps_st1[0][2*BIT_WIDTH_ML-1]=1;
/*
    for(int xt=0;xt< BIT_WIDTH_ML;xt++)
    {
        for(int yt=(2*BIT_WIDTH_ML)-1;yt>=0;yt--)
        {
                    printf("%d",pps_st1[xt][yt]);
        }
        printf("\n");
    }
*/

    int v_cnt=0;
    int s_cnt=0;
    int s_cnt1=0;

    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=0;
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]==1)
                {
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==2)
                {
                    pp_n[pp_cls]-=2;
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls] | pps_st1[v_cnt++][pp_cls];
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],&pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    if(c2==3)
                    {
                        approximate_compressor_4_2_c_3(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    else if(c2==4)
                    {
                        approximate_compressor_4_2_c_4(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n[pp_cls]=s_cnt;
    }

    pp_n[9]+=1;
    //stage 2
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l[pp_cls];
            s_cnt1=pp_n_l[pp_cls+1];
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]<=2)
                {
                    pps_st3[s_cnt++][pp_cls]=pps_st2[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],&pps_st3[s_cnt++][pp_cls], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], &pps_st3[s_cnt1++][pp_cls+1], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n_l[pp_cls]=s_cnt;
            pp_n_l[pp_cls+1]=s_cnt1;
    }

        //RCA on pps_st1 ->Series of full adders
    for(int pp_cls=0; pp_cls < ((2*BIT_WIDTH_ML)-k); pp_cls++)
    {
        pp_out=0;
        full_adder(pps_st3[0][pp_cls], pps_st3[1][pp_cls], cin, &pp_out,&cin);
        pp_product|=(long long)(pp_out<<pp_cls);
    }
    pp_product+=(long long)(pps_st3[0][((2*BIT_WIDTH_ML)-k)] << ((2*BIT_WIDTH_ML)-k));
    pp_product+=(long long)(cin << ((2*BIT_WIDTH_ML)-k));
    pp_product&=(1 << 2*BIT_WIDTH_ML)-1;
    pp_product=unsigned_to_signed(pp_product,2*BIT_WIDTH_ML);
    return(pp_product);
}


long long multilevel_opt(long long multiplicand, long long multiplier, char c2, char k_paper)
{

    char k=k_paper+2;
    uchar row_stop=0;
    long long pps_acc[BIT_WIDTH_ML]={ };
    long long mask=0;
    long long pps[BIT_WIDTH_ML]={ };
    long long pps_r=0;
    long long pp_product=0;
    long long pps_st2[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    long long pps_st3[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    long long pps_st1[BIT_WIDTH_ML][2*BIT_WIDTH_ML]={ };
    uchar carry[BIT_WIDTH_ML/4][2*BIT_WIDTH_ML]={ };
    uchar pp_n[2*BIT_WIDTH_ML]={ };
    uchar pp_n_l[2*BIT_WIDTH_ML]={ };
    int cnt=0;
    uchar pp_out=0;
    int pp_row1=0;
    uchar cin=0;
    long long multiplicand_in=multiplicand;
    long long multiplier_in=multiplier;
    //accurate part: NOTE it also requires the carry from the previous part
    char r_shift=0;
    for(int pp_row=0; pp_row<BIT_WIDTH_ML ;pp_row++)
    {
        pps_r=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        //complement the leading bits using XOR gate

        r_shift=(2*BIT_WIDTH_ML)-k-pp_row;
        if(r_shift<0)
        {
            pps_acc[pp_row]=pps_r<<((2*BIT_WIDTH_ML)-k-r_shift);
        }
        else{
                pps_acc[pp_row]=(pps_r>>r_shift)<<((2*BIT_WIDTH_ML)-k);
        }

        //print_binary2(pps_acc[pp_row],16);
        //printf("\n");

        pp_product+=pps_acc[pp_row];
    }

    //Have to shift carry fromt he approximate part by << (2*BIT_WIDTH_ML)-k;
    pp_product&=(1 << 2*BIT_WIDTH_ML)-1;
   // print_binary2(pp_product,16);
   // printf("\n");
   // return(pp_product);

    //approximate part
    for(int i=1;i<=BIT_WIDTH_ML;i++)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int i=BIT_WIDTH_ML-1;i>0;i--)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int pp_row=0; pp_row<BIT_WIDTH_ML ;pp_row++)
    {
        pps[pp_row]=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        pp_row1=pp_row;
        for(int pp_col=pp_row; pp_col < BIT_WIDTH_ML+pp_row;pp_col++)
        {
            if(pp_col>=BIT_WIDTH_ML)
            {
                pps_st1[--pp_row1][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
                //pp_row1--;
            }else
            {
                pps_st1[pp_row][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
            }

        }

    }

    //insert 1 at he BIT_WIDTH_ML+1 position ALSO increase the counter at this position.
/*
    for(int xt=0;xt< BIT_WIDTH_ML;xt++)
    {
        for(int yt=(2*BIT_WIDTH_ML)-1;yt>=0;yt--)
        {
                    printf("%d",pps_st1[xt][yt]);
        }
        printf("\n");
    }
*/

    int v_cnt=0;
    int s_cnt=0;
    int s_cnt1=0;

    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=0;
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]==1)
                {
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==2)
                {
                    pp_n[pp_cls]-=2;
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls] | pps_st1[v_cnt++][pp_cls];
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],&pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    if(c2==3)
                    {
                        approximate_compressor_4_2_c_3(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    else if(c2==4)
                    {
                        approximate_compressor_4_2_c_4(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                pp_n[pp_cls]-=4;
                }

            }
            pp_n[pp_cls]=s_cnt;
    }

    pp_n[9]+=1;
    //stage 2
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l[pp_cls];
            s_cnt1=pp_n_l[pp_cls+1];
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]<=2)
                {
                    pps_st3[s_cnt++][pp_cls]=pps_st2[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],&pps_st3[s_cnt++][pp_cls], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], &pps_st3[s_cnt1++][pp_cls+1], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n_l[pp_cls]=s_cnt;
            pp_n_l[pp_cls+1]=s_cnt1;
    }

        //RCA on pps_st1 ->Series of full adders
    for(int pp_cls=0; pp_cls < ((2*BIT_WIDTH_ML)-k); pp_cls++)
    {
        pp_out=0;
        full_adder(pps_st3[0][pp_cls], pps_st3[1][pp_cls], cin, &pp_out,&cin);
        pp_product|=(long long)(pp_out<<pp_cls);
    }
    pp_product+=(long long)(pps_st3[0][((2*BIT_WIDTH_ML)-k)] << ((2*BIT_WIDTH_ML)-k));
    pp_product+=(long long)(cin << ((2*BIT_WIDTH_ML)-k));
    pp_product&=(1 << 2*BIT_WIDTH_ML)-1;
    return(pp_product);
}

long long multilevel16_opt(long long multiplicand, long long multiplier, char c2, char k_paper)
{

    char k=k_paper+2;
    uchar row_stop=0;
    long long pps_acc[BIT_WIDTH_2ML]={ };
    long long mask=0;
    long long pps[BIT_WIDTH_2ML]={ };
    long long pps_r=0;
    long long pp_product=0;
    long long pps_st1[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st2[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st3[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st4[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    uchar carry[BIT_WIDTH_2ML/4][2*BIT_WIDTH_2ML]={ };
    uchar pp_n[2*BIT_WIDTH_2ML]={ };
    uchar pp_n_l[2*BIT_WIDTH_2ML]={ };
    uchar pp_n_l_2[2*BIT_WIDTH_2ML]={ };
    int cnt=0;
    uchar pp_out=0;
    int pp_row1=0;
    uchar cin=0;
    long long multiplicand_in=multiplicand;
    long long multiplier_in=multiplier;
    //accurate part: NOTE it also requires the carry from the previous part
    char r_shift=0;
    for(int pp_row=0; pp_row<BIT_WIDTH_2ML ;pp_row++)
    {
        pps_r=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        //complement the leading bits using XOR gate

        r_shift=(2*BIT_WIDTH_2ML)-k-pp_row;
        if(r_shift<0)
        {
            pps_acc[pp_row]=pps_r<<((2*BIT_WIDTH_2ML)-k-r_shift);
        }
        else{
                pps_acc[pp_row]=(pps_r>>r_shift)<<((2*BIT_WIDTH_2ML)-k);
        }

        //print_binary2(pps_acc[pp_row],16);
        //printf("\n");

        pp_product+=pps_acc[pp_row];
    }

    //Have to shift carry fromt he approximate part by << (2*BIT_WIDTH_2ML)-k;
    pp_product&=(1 << 2*BIT_WIDTH_2ML)-1;
   // print_binary2(pp_product,16);
   // printf("\n");
   // return(pp_product);

    //approximate part
    for(int i=1;i<=BIT_WIDTH_2ML;i++)
    {
        pp_n[cnt]=i;
        cnt++;

    }

    for(int i=BIT_WIDTH_2ML-1;i>0;i--)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int pp_row=0; pp_row<BIT_WIDTH_2ML ;pp_row++)
    {
        pps[pp_row]=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        pp_row1=pp_row;
        for(int pp_col=pp_row; pp_col < BIT_WIDTH_2ML+pp_row;pp_col++)
        {
            if(pp_col>=BIT_WIDTH_2ML)
            {
                pps_st1[--pp_row1][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
                //pp_row1--;
            }else
            {
                pps_st1[pp_row][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
            }

        }

    }

    //insert 1 at he BIT_WIDTH_2ML+1 position ALSO increase the counter at this position.
/*
    for(int xt=0;xt< BIT_WIDTH_2ML;xt++)
    {
        for(int yt=(2*BIT_WIDTH_2ML)-1;yt>=0;yt--)
        {
                    printf("%d",pps_st1[xt][yt]);
        }
        printf("\n");
    }
*/

    int v_cnt=0;
    int s_cnt=0;
    int s_cnt1=0;

    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=0;
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]==1)
                {
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==2)
                {
                    pp_n[pp_cls]-=2;
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls] | pps_st1[v_cnt++][pp_cls];
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],&pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    if(c2==3)
                    {
                        approximate_compressor_4_2_c_3(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    else if(c2==4)
                    {
                        approximate_compressor_4_2_c_4(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                pp_n[pp_cls]-=4;
                }

            }
            pp_n[pp_cls]=s_cnt;
    }

    pp_n[8]+=3;
    pp_n[9]+=2;
    pp_n[10]+=2;
    pp_n[11]+=2;

//stage 2
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l[pp_cls];
            s_cnt1=pp_n_l[pp_cls+1];
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]<=2)
                {
                    pps_st3[s_cnt++][pp_cls]=pps_st2[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],&pps_st3[s_cnt++][pp_cls], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], &pps_st3[s_cnt1++][pp_cls+1], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n_l[pp_cls]=s_cnt;
            pp_n_l[pp_cls+1]=s_cnt1;
    }

 //stage 2
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l_2[pp_cls];
            s_cnt1=pp_n_l_2[pp_cls+1];
            while(pp_n_l[pp_cls]>0)
            {
                if(pp_n_l[pp_cls]<=2)
                {
                    pps_st4[s_cnt++][pp_cls]=pps_st3[v_cnt++][pp_cls];
                    pp_n_l[pp_cls]-=1;
                }
                if(pp_n_l[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st3[v_cnt++][pp_cls],pps_st3[v_cnt++][pp_cls],pps_st3[v_cnt++][pp_cls],&pps_st4[s_cnt++][pp_cls], &pps_st4[s_cnt++][pp_cls]);
                    pp_n_l[pp_cls]-=3;
                }

                if(pp_n_l[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], &pps_st4[s_cnt1++][pp_cls+1], &pps_st4[s_cnt++][pp_cls]);
                    pp_n_l[pp_cls]-=4;
                }

            }
            pp_n_l_2[pp_cls]=s_cnt;
            pp_n_l_2[pp_cls+1]=s_cnt1;
    }

        //RCA on pps_st1 ->Series of full adders
    for(int pp_cls=0; pp_cls < ((2*BIT_WIDTH_2ML)-k); pp_cls++)
    {
        pp_out=0;
        full_adder(pps_st4[0][pp_cls], pps_st4[1][pp_cls], cin, &pp_out,&cin);
        pp_product|=(long long)(pp_out<<pp_cls);
    }
    pp_product+=(long long)(pps_st4[0][((2*BIT_WIDTH_2ML)-k)] << ((2*BIT_WIDTH_2ML)-k));
    pp_product+=(long long)(cin << ((2*BIT_WIDTH_2ML)-k));
    pp_product&=(1 << 2*BIT_WIDTH_2ML)-1;
    return(pp_product);
}


long long multilevel16_signed_opt(long long multiplicand_sign, long long multiplier_sign, char c2, char k_paper)
{
    long long factor=(long long)(pow(2,BIT_WIDTH_ML)-1);
    long long multiplicand=multiplicand_sign & factor;
    long long multiplier=multiplier_sign & factor;

    char k=k_paper+2;
    uchar row_stop=0;
    long long pps_acc[BIT_WIDTH_2ML]={ };
    long long mask=0;
    long long pps[BIT_WIDTH_2ML]={ };
    long long pps_r=0;
    long long pp_product=0;
    long long pps_st2[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st3[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st4[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st5[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    long long pps_st1[BIT_WIDTH_2ML][2*BIT_WIDTH_2ML]={ };
    uchar carry[BIT_WIDTH_2ML/4][2*BIT_WIDTH_2ML]={ };
    uchar pp_n[2*BIT_WIDTH_2ML]={ };
    uchar pp_n_l[2*BIT_WIDTH_2ML]={ };
    uchar pp_n_l_2[2*BIT_WIDTH_2ML]={ };
    uchar pp_n_l_3[2*BIT_WIDTH_2ML]={ };
    int cnt=0;
    uchar pp_out=0;
    int pp_row1=0;
    uchar cin=0;
    long long multiplicand_in=multiplicand;
    long long multiplier_in=multiplier;
    //accurate part: NOTE it also requires the carry from the previous part
    char r_shift=0;
    for(int pp_row=0; pp_row<BIT_WIDTH_2ML ;pp_row++)
    {
        pps_r=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        //complement the leading bits using XOR gate

        if(pp_row==BIT_WIDTH_2ML-1)
        {
            mask=(1 << (BIT_WIDTH_2ML-1));
            mask=(mask-1) |  (1 << (BIT_WIDTH_2ML));
            pps_r=pps_r ^ mask;
        }
        else if(pp_row==0)
        {
            pps_r=pps_r ^ (3 << (BIT_WIDTH_2ML-1));
        }
        else
        {
            pps_r=pps_r ^ (1 << (BIT_WIDTH_2ML-1));
        }
//        pps_acc[pp_row]=pps_r<<pp_row; //ONLY for test
        r_shift=(2*BIT_WIDTH_2ML)-k-pp_row;
        if(r_shift<0)
        {
            pps_acc[pp_row]=pps_r<<((2*BIT_WIDTH_2ML)-k-r_shift);
        }
        else{
                pps_acc[pp_row]=(pps_r>>r_shift)<<((2*BIT_WIDTH_2ML)-k);
        }

        //print_binary2(pps_acc[pp_row],16);
        //printf("\n");

        pp_product+=pps_acc[pp_row];
    }

    //Have to shift carry fromt he approximate part by << (2*BIT_WIDTH_2ML)-k;
    pp_product&=(1 << 2*BIT_WIDTH_2ML)-1;
   // print_binary2(pp_product,16);
   // printf("\n");
   // return(pp_product);

    //approximate part
    for(int i=1;i<=BIT_WIDTH_2ML;i++)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int i=BIT_WIDTH_2ML-1;i>0;i--)
    {
        pp_n[cnt]=i;
        cnt++;
    }

    for(int pp_row=0; pp_row<BIT_WIDTH_2ML ;pp_row++)
    {
        pps[pp_row]=((multiplier_in>>pp_row) & 1) * multiplicand_in;
        if(pp_row==BIT_WIDTH_2ML-1)
        {
            mask=(1 << (BIT_WIDTH_2ML-1));
            mask=(mask-1) |  (1 << (BIT_WIDTH_2ML));
            pps[pp_row]=pps[pp_row] ^ mask;
        }
        else if(pp_row==0)
        {
            pps[pp_row]=pps[pp_row] ^ (3 << (BIT_WIDTH_2ML-1));
        }
        else
        {
            pps[pp_row]=pps[pp_row] ^ (1 << (BIT_WIDTH_2ML-1));
        }
        pp_row1=pp_row;
        for(int pp_col=pp_row; pp_col < BIT_WIDTH_2ML+pp_row;pp_col++)
        {
            if(pp_col>=BIT_WIDTH_2ML)
            {
                pps_st1[--pp_row1][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
                //pp_row1--;
            }else
            {
                pps_st1[pp_row][pp_col]=(pps[pp_row] >> (pp_col-pp_row)) & 1;
            }

        }

    }

    //insert 1 at he BIT_WIDTH_2ML+1 position ALSO increase the counter at this position.
    pp_n[BIT_WIDTH_2ML]+=1;
    for(int pp_row=BIT_WIDTH_2ML-1;pp_row>0;pp_row--)
    {
        pps_st1[pp_row][BIT_WIDTH_2ML]=pps_st1[pp_row-1][BIT_WIDTH_2ML];
    }
    pps_st1[0][BIT_WIDTH_2ML]=1;
    pps_st1[0][2*BIT_WIDTH_2ML-1]=1;
/*
    for(int xt=0;xt< BIT_WIDTH_2ML;xt++)
    {
        for(int yt=(2*BIT_WIDTH_2ML)-1;yt>=0;yt--)
        {
                    printf("%d",pps_st1[xt][yt]);
        }
        printf("\n");
    }
*/

    int v_cnt=0;
    int s_cnt=0;
    int s_cnt1=0;

    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=0;
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]==1)
                {
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==2)
                {
                    pp_n[pp_cls]-=2;
                    pps_st2[s_cnt++][pp_cls]=pps_st1[v_cnt++][pp_cls] | pps_st1[v_cnt++][pp_cls];
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],pps_st1[v_cnt++][pp_cls],&pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    if(c2==3)
                    {
                        approximate_compressor_4_2_c_3(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    else if(c2==4)
                    {
                        approximate_compressor_4_2_c_4(pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], pps_st1[v_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls], &pps_st2[s_cnt++][pp_cls]);
                    }
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n[pp_cls]=s_cnt;
    }

//    pp_n[17]+=1;
  //  pp_n[18]+=1;
    //stage 2
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l[pp_cls];
            s_cnt1=pp_n_l[pp_cls+1];
            while(pp_n[pp_cls]>0)
            {
                if(pp_n[pp_cls]<=2)
                {
                    pps_st3[s_cnt++][pp_cls]=pps_st2[v_cnt++][pp_cls];
                    pp_n[pp_cls]-=1;
                }
                if(pp_n[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],pps_st2[v_cnt++][pp_cls],&pps_st3[s_cnt++][pp_cls], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=3;
                }

                if(pp_n[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], pps_st2[v_cnt++][pp_cls], &pps_st3[s_cnt1++][pp_cls+1], &pps_st3[s_cnt++][pp_cls]);
                    pp_n[pp_cls]-=4;
                }

            }
            pp_n_l[pp_cls]=s_cnt;
            pp_n_l[pp_cls+1]=s_cnt1;
    }

    //stage 3
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l_2[pp_cls];
            s_cnt1=pp_n_l_2[pp_cls+1];
            while(pp_n_l[pp_cls]>0)
            {
                if(pp_n_l[pp_cls]<=2)
                {
                    pps_st4[s_cnt++][pp_cls]=pps_st3[v_cnt++][pp_cls];
                    pp_n_l[pp_cls]-=1;
                }
                if(pp_n_l[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st3[v_cnt++][pp_cls],pps_st3[v_cnt++][pp_cls],pps_st3[v_cnt++][pp_cls],&pps_st4[s_cnt++][pp_cls], &pps_st4[s_cnt++][pp_cls]);
                    pp_n_l[pp_cls]-=3;
                }

                if(pp_n_l[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], pps_st3[v_cnt++][pp_cls], &pps_st4[s_cnt1++][pp_cls+1], &pps_st4[s_cnt++][pp_cls]);
                    pp_n_l[pp_cls]-=4;
                }

            }
            pp_n_l_2[pp_cls]=s_cnt;
            pp_n_l_2[pp_cls+1]=s_cnt1;
    }

   //stage 4
    for(int pp_cls=0;pp_cls < ((2*BIT_WIDTH_2ML)-k);pp_cls++)
    {
            v_cnt=0;
            s_cnt=pp_n_l_3[pp_cls];
            s_cnt1=pp_n_l_3[pp_cls+1];
            while(pp_n_l_2[pp_cls]>0)
            {
                if(pp_n_l_2[pp_cls]<=2)
                {
                    pps_st5[s_cnt++][pp_cls]=pps_st4[v_cnt++][pp_cls];
                    pp_n_l_2[pp_cls]-=1;
                }
                if(pp_n_l_2[pp_cls]==3)
                {
                    approximate_compressor_3_2(pps_st4[v_cnt++][pp_cls],pps_st4[v_cnt++][pp_cls],pps_st4[v_cnt++][pp_cls],&pps_st5[s_cnt++][pp_cls], &pps_st5[s_cnt++][pp_cls]);
                    pp_n_l_2[pp_cls]-=3;
                }

                if(pp_n_l_2[pp_cls]>=4)
                {
                    approximate_compressor_4_2_c_1(pps_st4[v_cnt++][pp_cls], pps_st4[v_cnt++][pp_cls], pps_st4[v_cnt++][pp_cls], pps_st4[v_cnt++][pp_cls], &pps_st5[s_cnt1++][pp_cls+1], &pps_st5[s_cnt++][pp_cls]);
                    pp_n_l_2[pp_cls]-=4;
                }

            }
            pp_n_l_3[pp_cls]=s_cnt;
            pp_n_l_3[pp_cls+1]=s_cnt1;
    }

        //RCA on pps_st1 ->Series of full adders
    for(int pp_cls=0; pp_cls < ((2*BIT_WIDTH_2ML)-k); pp_cls++)
    {
        pp_out=0;
        full_adder(pps_st5[0][pp_cls], pps_st5[1][pp_cls], cin, &pp_out,&cin);
        pp_product|=(long long)(pp_out<<pp_cls);
    }
    pp_product+=(long long)(pps_st5[0][((2*BIT_WIDTH_2ML)-k)] << ((2*BIT_WIDTH_2ML)-k));
    pp_product+=(long long)(cin << ((2*BIT_WIDTH_2ML)-k));
    pp_product&=(1 << 2*BIT_WIDTH_2ML)-1;
    pp_product=unsigned_to_signed(pp_product,2*BIT_WIDTH_2ML);
    return(pp_product);
}

//TESTED:OK
long long multilevel2B_opt(long long multiplicand, long long multiplier, char c2, char k_paper)
{
    long long multiplicand_l=0;
    long long multiplicand_h=0;
    long long multiplier_l=0;
    long long multiplier_h=0;

    long long prod1=0;
    long long prod2_l=0;
    long long prod2_h=0;
    long long prod2=0;
    long long prod3=0;
    long long prod=0;
    long long mask=0;

    mask= (1 << BIT_WIDTH_ML)- 1;
    multiplicand_l= multiplicand & mask;
    multiplicand_h= (multiplicand >> BIT_WIDTH_ML) & mask;
    multiplier_l= multiplier & mask;
    multiplier_h= (multiplier >> BIT_WIDTH_ML) & mask;

    prod1=multilevel_opt(multiplicand_h, multiplier_l,c2,k_paper);
    prod2_l=multilevel_opt(multiplicand_l, multiplier_l,c2,k_paper);
    prod2_h=multilevel_opt(multiplicand_h, multiplier_h,c2,k_paper);
    prod3=multilevel_opt(multiplicand_l, multiplier_h,c2,k_paper);

/*
    prod1=multiplicand_h*multiplier_l;
    prod2_l=multiplicand_l*multiplier_l;
    prod2_h=multiplicand_h* multiplier_h;
    prod3=multiplicand_l* multiplier_h;
*/
    prod1=prod1 << (BIT_WIDTH_ML);
    prod2=(prod2_h << (2*BIT_WIDTH_ML)) | prod2_l;
    prod3=prod3 << (BIT_WIDTH_ML);

    prod= prod1 + prod2 + prod3;

    return(prod);
}



//Test this by inserting a 8 bit accurate multiplier in place to approimate multiplier


long long multilevel2B_signed_opt(long long multiplicand_in, long long multiplier_in, char c2, char k_paper)
{
    long long multiplicand=0;
    long long multiplier=0;
    long long multiplicand_l=0;
    long long multiplicand_h=0;
    long long multiplier_l=0;
    long long multiplier_h=0;

    long long prod1=0;
    long long prod2_l=0;
    long long prod2_h=0;
    long long prod2=0;
    long long prod3=0;
    long long prod=0;
    long long mask=0;

    long long sign_bit=(multiplicand_in >> (2*BIT_WIDTH_ML)) ^ (multiplier_in >> (2*BIT_WIDTH_ML));
    sign_bit=sign_bit & 1;
    multiplicand=abs(multiplicand_in);
    multiplier=abs(multiplier_in);

    prod=multilevel2B_opt(multiplicand,multiplier,c2,k_paper) ;

    if(sign_bit==1)
    {
        prod=(~prod);
    }
    prod&=(1<<(4*BIT_WIDTH_ML))-1;
    prod|=sign_bit << (4*BIT_WIDTH_ML);
    prod=unsigned_to_signed(prod,4*BIT_WIDTH_ML+1);
    return(prod);
}

//AXBM
char encoder(uchar a, uchar b, uchar c, uchar d,uchar encoder_type)
{

    char factor=(-4*a)+(2*b)+(c)+(d);

    if(encoder_type==LUT_AXBM_EN_1)
    {

        if(a==0 && b==0 && c==0 && d==0)
        {
            factor=1;
        }
        if(a==1 && b==1 && c==1 && d==1)
        {
            factor=-1;
        }

        if(a==0 && b==1 && c==0 && d==1)
        {
            factor=4;
        }
        if(a==0 && b==1 && c==1 && d==0)
        {
            factor=4;
        }

        if(a==1 && b==0 && c==0 && d==1)
        {
            factor=-4;
        }
        if(a==1 && b==0 && c==1 && d==0)
        {
            factor=-4;
        }
    }

    if(encoder_type==LUT_AXBM_EN_2)
    {

        if(a==0 && b==0 && c==0 && d==0)
        {
            factor=1;
        }
        if(a==1 && b==1 && c==1 && d==1)
        {
            factor=-1;
        }

        if(a==0 && b==1 && c==0 && d==1)
        {
            factor=4;
        }
        if(a==0 && b==1 && c==1 && d==0)
        {
            factor=4;
        }

        if(a==1 && b==0 && c==0 && d==1)
        {
            factor=-2;
        }
        if(a==1 && b==0 && c==1 && d==0)
        {
            factor=-2;
        }
    }


    return(factor);

}

long long ax_bm1_opt(long multiplicand, long multiplier)
{
    //long long multiplier= multiplier_in & ((1 << BIT_WIDTH)-1);
    long multiplier_appended=multiplier<<1;
    long multiplicand_appended=multiplicand<<2;
    uchar bits[4]={ };
    long long pp_products[PP_NUM_UNSIGNED_AXBM]={ };
    long long prod=0;
    int lut_counter=0;
    uchar multiple_tuple[3];
    long long e=0;
    long long header=0;
    //condition t_row
    uchar encode[2]={ };
    uchar pp_result[2]={ };

    for(int pp_row=0;pp_row<PP_NUM_UNSIGNED_AXBM;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(3*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(3*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(3*pp_row+2) & 1;
        bits[BIT_MSB1]=multiplier_appended >>(3*pp_row+3) & 1;

        pp_products[pp_row]=encoder(bits[BIT_MSB1],bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],LUT_AXBM_EN_1) * multiplicand;
    }

    for(int pp_row=0;pp_row<PP_NUM_UNSIGNED_AXBM;pp_row++)
    {
       prod+=(pp_products[pp_row])<<(3*pp_row);
    }

    return(prod);
}

long long ax_bm2_opt(long multiplicand, long multiplier)
{
    //long long multiplier= multiplier_in & ((1 << BIT_WIDTH)-1);
    long multiplier_appended=multiplier<<1;
    long multiplicand_appended=multiplicand<<2;
    uchar bits[4]={ };
    long long pp_products[PP_NUM_UNSIGNED_AXBM]={ };
    long long prod=0;
    int lut_counter=0;
    uchar multiple_tuple[3];
    long long e=0;
    long long header=0;
    //condition t_row
    uchar encode[2]={ };
    uchar pp_result[2]={ };

    for(int pp_row=0;pp_row<PP_NUM_UNSIGNED_AXBM;pp_row++)
    {
        //2i-1, 2i, 2i+1; However since the index of the multiplier_appended starts from 0 instead of -1, just +1 for all indices
        bits[BIT_LSB]=multiplier_appended >>(3*pp_row) & 1;
        bits[BIT_CSB]=multiplier_appended >>(3*pp_row+1) & 1;
        bits[BIT_MSB]=multiplier_appended >>(3*pp_row+2) & 1;
        bits[BIT_MSB1]=multiplier_appended >>(3*pp_row+3) & 1;

        pp_products[pp_row]=encoder(bits[BIT_MSB1],bits[BIT_MSB],bits[BIT_CSB],bits[BIT_LSB],LUT_AXBM_EN_2) * multiplicand;
    }

    for(int pp_row=0;pp_row<PP_NUM_UNSIGNED_AXBM;pp_row++)
    {
       prod+=(pp_products[pp_row])<<(3*pp_row);
    }
    //post pocessing to truncate 9 LSBs in AX_BM2
    if(BIT_WIDTH==16)
    {
        prod&=(-1<<9);
        prod=prod+(long long)(1<<9);
    }
    return(prod);
}

//Funtion to be called in python
long long approximate_multiplier(long long multiplicand, long long multiplier, uchar MULT)
{
    long long result=0;
    switch(MULT)
    {
        case BOOTH_APPROX:
            result=area_optimized_hdl_approx(multiplicand,multiplier);
            break;
        case MUL8_13_4:
            result=multilevel_signed_opt(multiplicand, multiplier,3,4);
            break;
        case MUL8_14_4:
            result=multilevel_signed_opt(multiplicand, multiplier,4,4);
            break;
        case LC_BW_2:
            result=lc_baugh_wooley_opt(multiplicand,multiplier,2);
            break;
        case LC_BW_2_AC:
            result=lc_baugh_wooley_opt_ac(multiplicand,multiplier,2);
            break;
        case LC_BOOTH_2:
            result=lc_booth_opt(multiplicand,multiplier,2);
            break;
        case LC_BOOTH_2_AC:
            result=lc_booth_opt_ac(multiplicand,multiplier,2);
            break;
        case AXBM1:
            result=ax_bm1_opt(multiplicand,multiplier);
            break;
        case AXBM2:
            result=ax_bm2_opt(multiplicand,multiplier);
            break;
        default:
            printf("Invalid choice\n");
            break;

    }
    return(result);
}

/*
long long approximate16_multiplier(long long multiplicand, long long multiplier, uchar MULT)
{
    long long result=0;
    switch(MULT)
    {
        case BOOTH_APPROX:
            result=area_optimized_hdl_approx(multiplicand,multiplier);
            break;
        case MUL8_13_4:
            result=multilevel2B_signed3_opt(multiplicand, multiplier,3,4);
            break;
        case MUL8_14_4:
            result=multilevel2B_signed3_opt(multiplicand, multiplier,4,4);
            break;
        case LC_BW_2:
            result=lc_baugh_wooley_opt(multiplicand,multiplier,2);
            break;
        case LC_BW_2_AC:
            result=lc_baugh_wooley_opt_ac(multiplicand,multiplier,2);
            break;
        case LC_BOOTH_2:
            result=lc_booth_opt(multiplicand,multiplier,2);
            break;
        case LC_BOOTH_2_AC:
            result=lc_booth_opt_ac(multiplicand,multiplier,2);
            break;
        case AXBM1:
            result=ax_bm1_opt(multiplicand,multiplier);
            break;
        case AXBM2:
            result=ax_bm2_opt(multiplicand,multiplier);
            break;
        default:
            printf("Invalid choice\n");
            break;

    }
    return(result);
}*/

double ber_error(long long acc, long long app,char p_bit_width)
{
        double error;
        char sum=0;
        long long error_stream= acc ^ app;
        for(int i=0;i<p_bit_width;i++)
        {
            sum+=(error_stream >> i) & 1;
        }
        error=(double)sum/(double)p_bit_width;
        return(error);

}

void error_benchmarks(double* errors)
{
    //get absolute value of pa-pt
    //operate on this value to get error metrics.
    //change or make 'unsigned int' generic since it holds result of Booth mutliplication.
    // SHould a * b  and b* a considered separately??? SHould it be C^n_2 (unique comninations)? DOesn't matter since 2*(same_error) /2 =same_error
    long long M_acc=0;
    long long M_app=0;
    long long ed=0;
    double red=0;
    double redv=0;
    double max_re=0;
    long long range=0;
    long long mred_range=0;
    long long error_counter=0;
    long long max_error=0;
    long long M_acc_max_value=0; //should be signed maximum or abs maximum???
    double complete_percentage=0;
    double ber_sum=0;
    double bit_probabilities[2*BIT_WIDTH]={ };
    double ln_prob_0[2*BIT_WIDTH]={ };
    double ln_prob_1[2*BIT_WIDTH]={ };
    double pr=0;
    //double pr_num[1<<(2*BIT_WIDTH)]={ };
    long long number_elements=((long long)1)<<(2*BIT_WIDTH);
    long long products_num[((long long)1)<<(2*BIT_WIDTH)]={ };
    long cnt=0;
    for(long long a=SIGNED_LOWER_LIMIT; a<=SIGNED_UPPER_LIMIT;a++)
     {
        for(long long b=SIGNED_LOWER_LIMIT; b<=SIGNED_UPPER_LIMIT;b++)
        {
            // Accurate booth value = M_acc   ? Isn it just Multiplier x Multipicand???????????
            // Approximate both value =M_app
            //M_app=area_optimized_hdl_approx(a,b);
            //M_app=star_booth_opt_v3(a,b);
            //M_app=lc_booth_opt2(a,b,2);
//            M_app=lc_baugh_wooley_opt2(a,b,2);
            //M_app=unsigned_to_signed(M_app,BIT_WIDTH*2);
            //M_app=ax_bm1(a,b);
            //M_app=baugh_wooley_hdl(a,b);
            //M_app=baugh_wooley_hdl_v2(a,b,2);
            //M_app=multilevel_signed_opt(a,b,4,4);
            //M_app=unsigned_to_signed(M_app,BIT_WIDTH*2);
            //M_app=ax_bm1_opt(a,b);
            //M_app=unsigned_to_signed(M_app,BIT_WIDTH*2);
            //M_app=multilevel2B_opt(a,b,3,4);
            //M_app=star_booth_opt_v3(a,b);
          //  M_app=booth_opt2(a,b);
            //M_app=unsigned_to_signed(M_app,2*BIT_WIDTH);
            //M_app=approximate_multiplier(a,b,AXBM2);
            //M_app=multilevel2B_signed_opt(a,b,3,4);
            //M_app=ax_bm2_opt(a,b);
            //M_app=multilevel_opt(a,b,4,4);
            //M_app=multilevel2B_opt(a,b,4,4);
            //M_app=multilevel_signed_opt(a,b,3,4);
            //M_app=area_optimized_hdl_approx(a,b);
            //M_app=area_optimized_hdl_approx(a,b);
            //M_app=ax_bm2_opt(a,b);
            //M_app=lc_baugh_wooley_opt_ac(a,b,2);
            //M_app=lc_(a,b,2);
            //M_app=lc_booth_opt(a,b,4);
            //M_app=lc_baugh_wooley_opt(a,b,1);
            //products_num[cnt++]=(double)M_app;
            M_acc=a * b;
            products_num[cnt++]=(long long)M_acc;
            //compute_probabilities(M_app,bit_probabilities);
            ed=ed+abs(M_app-M_acc);
            //Push((double)abs(M_acc-M_app));
            ber_sum+=ber_error(M_acc,M_app,2*BIT_WIDTH);
            if(M_acc)
            {

                red=red+abs(M_app-M_acc)/(double)abs(M_acc);;
                if(max_re < (abs(M_app-M_acc)/(double)abs(M_acc))){
                    max_re=(abs(M_app-M_acc)/(double)abs(M_acc));
                }
                mred_range+=1;
            }
            if(max_error < abs(M_app-M_acc)){
                max_error=abs(M_app-M_acc);}

            if(M_acc_max_value < abs(M_acc)){
                M_acc_max_value=abs(M_acc);}

            if(M_app!=M_acc)
            {
                error_counter+=1;
            }
            range+=1;

        }
        complete_percentage=range*100;
        printf("%Lf\n",complete_percentage/(long double)(pow(2,2*BIT_WIDTH)));
    }
    double mean=Mean();
    double variance=Variance();
    double std=StandardDeviation();
    double med=ed/(range);
    errors[MED]=med;
    double mred=red/(mred_range);
    errors[MRED]=mred;
    double percentage_error=100*(error_counter)/range;
    //test
    //M_acc_max_value=pow((pow(2,BIT_WIDTH)-1),2);
    double nmed=med/(M_acc_max_value);
    //double nmed1=med/pow((pow(2,BIT_WIDTH)-1),2);
    double ber_error=ber_sum/(double)range;
    errors[NMED]=nmed;
    noramlize_probabilities((double)range,bit_probabilities);
    log_prob(ln_prob_0,ln_prob_1,bit_probabilities);
    //pr=compute_probablity_num(1000,ln_prob_0,ln_prob_1);
    //compute_probability_distribution(pr_num,ln_prob_0,ln_prob_1);
    //write_csv(pr_num,1<<(2*BIT_WIDTH),"one_app_4.csv");
    write_csv(products_num,number_elements,"multiplication_16_c.csv");
}




void test_encoder()
{
    int lut_counter=0;
    uchar encode[2]={ };
    for(uchar a=0;a<=1;a++)
    {
        for(uchar b=0;b<=1;b++)
        {
            for(uchar c=0;c<=1;c++)
            {
                for(uchar d=0;d<=1;d++)
                {
                        approximate_compressor_4_2_c_4(a,b,c,d,&encode[0],&encode[1]);
                        printf("a=%d b=%d c=%d d=%d  x1=%d x2=%d\n",a,b,c,d,encode[0],encode[1]);
                }
            }
            }
            }

}

//void generate_hex(uchar a, uchar b, uchar c, uchar d, uchar e, uchar f, uchar *o6, uchar mode)
long long generate_hex(uchar mode)
{
    char ptr=0;
    unsigned long long hex_val=0;
    uchar res[2]={ };
    int lut_counter=0;
    if(mode==1)
    {
        for(uchar a=0;a<2;a++)
            for(uchar b=0;b<2;b++)
                for(uchar c=0;c<2;c++)
                    for(uchar d=0;d<2;d++)
                        for(uchar e=0;e<2;e++)
                            for(uchar f=0;f<2;f++)
                            {
                                printf("a=%d b=%d c=%d d=%d e=%d f=%d\n",a,b,c,d,e,f);
                                LUT_6_2(a, b, c, d, e, f, LUT_TEST, &lut_counter, res);
                                hex_val=hex_val | (unsigned long long)(res[0]<<ptr);
                                ptr++;
                            }
    }
    return(hex_val);
}


void main()
{
    char s1=0;
    uchar c1=2;
    long long prod=0;
    unsigned long long hx_val=0;
    uchar carry1[PP_NUM_SIGNED*(BIT_WIDTH+2)];
    double err[3];
    error_benchmarks(&err);
    printf("NMED=%d",err[NMED]);
}
