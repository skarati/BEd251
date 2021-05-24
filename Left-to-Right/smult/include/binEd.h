#ifndef BINKL_H_
#define BINKL_H_
#include "binvecmult.h"

inline void clamp(unsigned char n[32]);
inline void invert(gfe1x *op, gfe1x *in);
extern void conditionalSwap(gfe1x *a, gfe1x *b, gfe1x *c, gfe1x *d, vec sbit);
inline void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]);
inline void scalar_mult_var_base(unsigned char op[32], unsigned char n[32], unsigned char varbase[32]);

inline void clamp(unsigned char n[32]){
	n[0] 	= n[0] & 0xfc;
	n[31] 	= n[31] & 0x07;
	n[31] 	= n[31] | 0x04;
}


void scalar_mult_fixed_base(unsigned char op[32], unsigned char n[32]){
	int 	prevbit, bit, i, j, k;
	vec 	vbit, t1;
	gfe1x 	sw,sz,rw,rz,vbx;
	gfe1x 	t2,t3,t4;
	gfe 	op_gfe;
	gfe	work;

	vbx.v[0] = baseP[0];   vbx.v[1] = baseP[1];   vbx.v[2] = baseP[2];   vbx.v[3] = baseP[3];  
	sw.v[0] = baseP[0];    sw.v[1] = baseP[1];    sw.v[2] = baseP[2];    sw.v[3] = baseP[3];
	sz.v[0] = one;         sz.v[1] = zero;        sz.v[2] = zero;        sz.v[3] = zero;
	rw.v[0] = _2baseP[0];  rw.v[1] = _2baseP[1];  rw.v[2] = _2baseP[2];  rw.v[3] = _2baseP[3];
	rz.v[0] = one;         rz.v[1] = zero;        rz.v[2] = zero;        rz.v[3] = zero;
	

	j = 2;	
	i = 31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1){
			j = 7;
			i--;
		}
	}


	prevbit = 0;
	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = ((n[i]>>j) & 1); 

			//Conditional Swap
			vbit = _mm_set1_epi32(0-(bit^prevbit));
			conditionalSwap(&sw, &rw, &sz, &rz, vbit);

			ladderStep(&sw, &sz, &rw, &rz, &vbx);
			prevbit = bit;
		}

		j=7;
	}


	//Conditional Swap
	vbit = _mm_set1_epi32(0-bit);
	conditionalSwap(&sw, &rw, &sz, &rz, vbit);

	invert(&sz,&sz);
	gfe1xMult(&sw, &sw, &sz);

	convert_gfe1x2gfe(&op_gfe, &sw);
	convert_itoc(&op_gfe, op);
}


void scalar_mult_var_base(unsigned char op[32], unsigned char n[32], unsigned char varbase[32]){
	int 	prevbit, bit, i, j, k;
	vec 	vbit, t1;
	gfe1x 	sw,sz,rw,rz,vbx;
	gfe1x 	t2,t3,t4;
	gfe 	op_gfe;
	gfe	work;

	convert_ctoi(&work,varbase);
	convert_gfe2gfe1x(&vbx,&work);

	sw = vbx;
	sz.v[0] = one; sz.v[1]=zero; sz.v[2]=zero; sz.v[3]=zero;
	
	//Set r = KL_Point_Double(s)
	gfe1xAdd(&t2, &sw, &sz);
	gfe1xMult(&t2, &sw, &t2);
	gfe1xSq(&rw,&t2);
	rz.v[0] = _mm_xor_si128(rw.v[0],d1); 
	rz.v[1]=rw.v[1]; rz.v[2]=rw.v[2]; rz.v[3]=rw.v[3];
				
	
	j = 2;	
	i = 31;
	while(bit == 0){
		bit = (n[i]>>j) & 1;
		j--;
		if(j==-1){
			j = 7;
			i--;
		}
	}


	prevbit = 0;
	for(;i>=0;i--){
    		for(;j>=0;j--){
			bit = ((n[i]>>j) & 1); 

			//Conditional Swap
			vbit = _mm_set1_epi32(0-(bit^prevbit));
			conditionalSwap(&sw, &rw, &sz, &rz, vbit);

			ladderStep(&sw, &sz, &rw, &rz, &vbx);
			prevbit = bit;
		}

		j=7;
	}


	//Conditional Swap
	vbit = _mm_set1_epi32(0-bit);
	conditionalSwap(&sw, &rw, &sz, &rz, vbit);

	invert(&sz,&sz);
	gfe1xMult(&sw, &sw, &sz);

	convert_gfe1x2gfe(&op_gfe, &sw);
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
