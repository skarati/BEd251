#ifndef BINKL_H_
#define BINKL_H_
#include "binvecmult.h"
#include "precomp.h"

inline void clamp(unsigned char n[32]);
inline void invert(gfe1x *op, gfe1x *in);
extern void conditionalSwap(gfe1x *a, gfe1x *b, gfe1x *c, gfe1x *d, vec sbit);
inline void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]);

inline void clamp(unsigned char n[32]){
	n[0] 	= n[0] & 0xfc;
	n[0]	= n[0] | 0x04;
	n[31] 	= n[31] & 0x07;
	n[31] 	= n[31] | 0x04;
}


align void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]){
	int 	prevbit, bit, i, j, k;
	vec 	vbit, t1;
	gfe1x 	r0w,r1w,r1z,r2w,r2z;
	gfe1x 	t2,t3,t4,t22;
	gfe 	op_gfe;

	

	r0w.v[0] = precompBase[1][0]; r0w.v[1] = precompBase[1][1]; r0w.v[2] = precompBase[1][2]; r0w.v[3] = precompBase[1][3];
	
	r1w.v[0] = zero;   r1w.v[1]=zero; r1w.v[2]=zero; r1w.v[3]=zero;
	r1z.v[0] = one;    r1z.v[1]=zero; r1z.v[2]=zero; r1z.v[3]=zero;
   
	r2w.v[0] = precompBase[1][0]; r2w.v[1] = precompBase[1][1]; r2w.v[2] = precompBase[1][2]; r2w.v[3] = precompBase[1][3];
	r2z.v[0] = one; r2z.v[1]=zero; r2z.v[2]=zero; r2z.v[3]=zero;


        k = 2;
	j=2;	
	prevbit = 1;
	for(i=0;i<31;i++){
    		for(;j<8;j++){
			bit = ((n[i]>>j) & 1); 
			
			//Conditional Swap
			vbit = _mm_set1_epi32(0-(bit^prevbit));
			conditionalSwap(&r1w, &r2w, &r1z, &r2z, vbit);
			
			ladderStep(&r0w, &r1w, &r1z, &r2w, &r2z);

        	        r0w.v[0] = precompBase[k][0]; r0w.v[1] = precompBase[k][1]; r0w.v[2] = precompBase[k][2]; r0w.v[3] = precompBase[k][3];
        	        k++;
			
			prevbit = bit;
		}
		j=0;
	}
	
	for(j=0;j<3;j++){
		bit = ((n[31]>>j) & 1); 

		//Conditional Swap
		vbit = _mm_set1_epi32(0-(bit^prevbit));
		conditionalSwap(&r1w, &r2w, &r1z, &r2z, vbit);
			
		ladderStep(&r0w, &r1w, &r1z, &r2w, &r2z);
                r0w.v[0] = precompBase[k][0]; r0w.v[1] = precompBase[k][1]; r0w.v[2] = precompBase[k][2]; r0w.v[3] = precompBase[k][3];
                k++;
		prevbit = bit;
	}


	invert(&r1z,&r1z);
	gfe1xMult(&r1w, &r1w, &r1z);

	convert_gfe1x2gfe(&op_gfe, &r1w);
	convert_itoc(&op_gfe, op);
}

inline void invert(gfe1x *op, gfe1x *in){
	gfe1x t;
	gfe1x x2,x3,x4,x7,x_6_1,x_12_1,x_24_1,x_25_1,x_50_1;
	gfe1x x_100_1,x_125_1,x_250_1;

	/*2*/			gfe1xSq(&x2, in);				//1S
	/*3*/			gfe1xMult(&x3, &x2, in);			//1M
	/*4*/			gfe1xSq(&x4, &x2);				//2S
	/*7*/			gfe1xMult(&x7, &x4, &x3);			//2M


	/*2^6-3*/		gfe1xnSq(&x_6_1, &x7,3);			//5S
	/*2^6-1*/		gfe1xMult(&x_6_1, &x_6_1,&x7);			//3M
	
	/*2^12-6*/		gfe1xnSq(&x_12_1, &x_6_1,6);			//11S
	/*2^12-1*/		gfe1xMult(&x_12_1, &x_12_1,&x_6_1);		//4M

	/*2^24-1*/		gfe1xnSq(&x_24_1, &x_12_1,12);			//23S	
	/*2^24-1*/		gfe1xMult(&x_24_1, &x_24_1,&x_12_1);		//5M
	
	/*2^25-1*/		gfe1xSq(&x_25_1, &x_24_1);			//24S	
	/*2^25-1*/		gfe1xMult(&x_25_1, &x_25_1,in);			//6M
	
	/*2^50-1*/		gfe1xnSq(&x_50_1, &x_25_1,25);			//49S
	/*2^50-1*/		gfe1xMult(&x_50_1,&x_50_1, &x_25_1);		//7M
	
	/*2^100-1*/		gfe1xnSq(&x_100_1, &x_50_1,50);			//99S
	/*2^100-1*/		gfe1xMult(&x_100_1,&x_100_1, &x_50_1);		//8M
	
	/*2^125-1*/		gfe1xnSq(&x_125_1, &x_100_1,25);		//124S
	/*2^125-1*/		gfe1xMult(&x_125_1,&x_125_1, &x_25_1);		//9M
	
	/*2^250-1*/		gfe1xnSq(&x_250_1, &x_125_1,125);		//249S
	/*2^250-1*/		gfe1xMult(&x_250_1,&x_250_1, &x_125_1);		//10M

	/*2^251-2*/		gfe1xSq(op,&x_250_1);				//250S

}

#endif
