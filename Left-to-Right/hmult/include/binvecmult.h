#ifndef BINVECMUL_H_
#define BINVECMUL_H_

#define vADD(C,A,B)    {C = _mm_xor_si128(A,B);}
#define vSUB(C,A,B)    vADD(C,A,B)
#define vMULT(C,A,B)   {C = _mm_clmulepi64_si128(A,B,0);}
#define vSFTLB(C,A,B)   {C= _mm_slli_si128(A,B);}
#define vSFTRB(C,A,B)   {C= _mm_srli_si128(A,B);}
#define vAND(C,A,B)    {C= _mm_and_si128(A,B);}


extern void bincopy(gfe1x *c, gfe1x *a);
extern void gfe1xnSq(gfe1x *c, gfe1x *a, int n);
extern void gfe1xMult(gfe1x *c, gfe1x *a, gfe1x *b);
extern void gfe1xMultConst(gfe1x *c, gfe1x *a, vec b);
extern void gfe1xAdd(gfe1x *c, gfe1x *a, gfe1x *b);
extern void ladderStep(gfe1x *sw, gfe1x *sz, gfe1x *rw, gfe1x *rz, gfe1x *vbx);
#define gfe1xSq(x,y) gfe1xnSq(x,y,1)


#endif

