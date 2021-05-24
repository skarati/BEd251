.p2align 5
.globl ladderStep
ladderStep:

	movq   %rsp, %r11			#r11 = rsp
	subq   $512, %rsp			#create a stack of 512 byte
	movq   %r11,  0(%rsp)			#push r11 in to stack to restore rsp at the end

	
	###############
	# A= w0*w1#
	###############
	movdqa	0(%rdi),%xmm8
	movdqa	16(%rdi),%xmm9
	movdqa	32(%rdi),%xmm10	
	movdqa	48(%rdi),%xmm11	

	movdqa	0(%rsi),%xmm12
	movdqa	16(%rsi),%xmm13	
	movdqa	32(%rsi),%xmm14	
	movdqa	48(%rsi),%xmm15	

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9
	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,8(%rsp)			#store A to stack
	movdqa	%xmm1,24(%rsp)			#store A to stack
	movdqa	%xmm2,40(%rsp)			#store A to stack
	movdqa	%xmm3,56(%rsp)			#store A to stack
	
	
	###############
	# B = w1+z1#
	###############
	movdqa	0(%rdx),%xmm8
	movdqa	16(%rdx),%xmm9
	movdqa	32(%rdx),%xmm10	
	movdqa	48(%rdx),%xmm11	

	movdqa	0(%rsi),%xmm12
	movdqa	16(%rsi),%xmm13	
	movdqa	32(%rsi),%xmm14	
	movdqa	48(%rsi),%xmm15	

	pxor %xmm8, %xmm12			#xmm12 = xmm8 ^ xmm12 = xmm8 + xmm12
	pxor %xmm9, %xmm13			#xmm13 = xmm9 ^ xmm13 = xmm9 + xmm13
	pxor %xmm10, %xmm14			#xmm14 = xmm10 ^ xmm14 = xmm10 + xmm14
	pxor %xmm11, %xmm15			#xmm15 = xmm11 ^ xmm15 = xmm11 + xmm15

	####################
	# B = (w0+1)(w1+z1)#
	####################
	movdqa	0(%rdi),%xmm8
	movdqa	16(%rdi),%xmm9
	movdqa	32(%rdi),%xmm10	
	movdqa	48(%rdi),%xmm11	

	pxor	one, %xmm8
	
	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9
	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,72(%rsp)			#store B to stack
	movdqa	%xmm1,88(%rsp)			#store B to stack
	movdqa	%xmm2,104(%rsp)			#store B to stack
	movdqa	%xmm3,120(%rsp)			#store B to stack	

	###############
	#t1 = A+B     #
	###############	

	movdqa	8(%rsp),%xmm12			#restore A to stack
	movdqa	24(%rsp),%xmm13			#restore A to stack
	movdqa	40(%rsp),%xmm14 		#restore A to stack
	movdqa	56(%rsp),%xmm15			#restore A to stack

	pxor %xmm0,%xmm12
	pxor %xmm1,%xmm13
	pxor %xmm2,%xmm14
	pxor %xmm3,%xmm15


	###############
	#t1 = A+B+z1  #
	###############	

	movdqa	0(%rdx),%xmm0		
	movdqa	16(%rdx),%xmm2		
	movdqa	32(%rdx),%xmm4 	
	movdqa	48(%rdx),%xmm6		

	pxor %xmm12,%xmm0
	pxor %xmm13,%xmm2
	pxor %xmm14,%xmm4
	pxor %xmm15,%xmm6


	#############
	#t1 = t1^2  #
	#############
	
	#movdqa %xmm3, %xmm6			#xmm6 = xmm3
	#movdqa %xmm2, %xmm4			#xmm4 = xmm2
	#movdqa %xmm1, %xmm2			#xmm2 = xmm1

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7
	
	
	#############
	#t1 = d1*t1 #
	#############
		
	pclmulqdq $0x00, d1, %xmm0		#xmm0 = xmm0 * d1
	pclmulqdq $0x00, d1, %xmm1		#xmm1 = xmm1 * d1
	pclmulqdq $0x00, d1, %xmm2		#xmm2 = xmm2 * d1
	pclmulqdq $0x00, d1, %xmm3		#xmm3 = xmm3 * d1

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7


	#############
	#w5 = z2*t1 #
	#############
	movdqa	%xmm0,%xmm8
	movdqa	%xmm1,%xmm9
	movdqa	%xmm2,%xmm10	
	movdqa	%xmm3,%xmm11	

	movdqa	0(%r8),%xmm12
	movdqa	16(%r8),%xmm13
	movdqa	32(%r8),%xmm14	
	movdqa	48(%r8),%xmm15	

	#movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9
	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,0(%rsi)			
	movdqa	%xmm1,16(%rsi)			
	movdqa	%xmm2,32(%rsi)			
	movdqa	%xmm3,48(%rsi)			


	#############
	#t2 = z1^2 #
	#############
	movdqa	0(%rdx),%xmm0		
	movdqa	16(%rdx),%xmm2		
	movdqa	32(%rdx),%xmm4 	
	movdqa	48(%rdx),%xmm6	

	pclmulqdq $0x00, %xmm0, %xmm0 		#Binary field multiplication based on imm8=$0x00 as xmm0 = xmm0*xmm0, output is stored at last xmm
	pclmulqdq $0x00, %xmm2, %xmm2		#computation of xmm2^2
	pclmulqdq $0x00, %xmm4, %xmm4		#computation of xmm4^2
	pclmulqdq $0x00, %xmm6, %xmm6		#computation of xmm6^2

	movdqa %xmm0, %xmm1			#set xmm1 = xmm0  
	movdqa %xmm2, %xmm3			#set xmm3 = xmm2  
	movdqa %xmm4, %xmm5			#set xmm5 = xmm4
	movdqa %xmm6, %xmm7			#set xmm7 = xmm6

	#expand
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm1			#xmm1 = xmm1>>64
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm3			#xmm3 = xmm3>>64
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5>>64
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7


	#############
	#t2 = d1*t2 #
	#############
		
	pclmulqdq $0x00, d1, %xmm0		#xmm0 = xmm0 * d1
	pclmulqdq $0x00, d1, %xmm1		#xmm1 = xmm1 * d1
	pclmulqdq $0x00, d1, %xmm2		#xmm2 = xmm2 * d1
	pclmulqdq $0x00, d1, %xmm3		#xmm3 = xmm3 * d1

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7


	movdqa	%xmm0,136(%rsp)			
	movdqa	%xmm1,152(%rsp)			
	movdqa	%xmm2,168(%rsp)			
	movdqa	%xmm3,184(%rsp)			
	
	
	#############
	#t3 = A*B   #
	#############
	movdqa	8(%rsp),%xmm8			#restore A from stack
	movdqa	24(%rsp),%xmm9			#restore A from stack
	movdqa	40(%rsp),%xmm10			#restore A from stack
	movdqa	56(%rsp),%xmm11			#restore A from stack
	
	movdqa	72(%rsp),%xmm12			#restore B from stack
	movdqa	88(%rsp),%xmm13			#restore B from stack
	movdqa	104(%rsp),%xmm14		#restore B from stack
	movdqa	120(%rsp),%xmm15		#restore B from stack

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9
	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	#############
	#t3 = t2+t3 #
	#############

	movdqa	136(%rsp),%xmm8			#restore t2 from stack
	movdqa	152(%rsp),%xmm9			#restore t2 from stack
	movdqa	168(%rsp),%xmm10		#restore t2 from stack
	movdqa	184(%rsp),%xmm11		#restore t2 from stack

	pxor %xmm0,%xmm8
	pxor %xmm1,%xmm9
	pxor %xmm2,%xmm10
	pxor %xmm3,%xmm11

	#############
	#z5 = x2*t3 #
	#############
	movdqa	0(%rcx),%xmm12			
	movdqa	16(%rcx),%xmm13			
	movdqa	32(%rcx),%xmm14			
	movdqa	48(%rcx),%xmm15			

	movdqa	%xmm8, %xmm0			#xmm0 = xmm8	
	pclmulqdq $0x00, %xmm12, %xmm0		#xmm0 = xmm8 * xmm12
	movdqa	%xmm9, %xmm2			#xmm2 = xmm9
	pclmulqdq $0x00, %xmm13, %xmm2		#xmm2 = xmm9 * xmm13
	movdqa 	%xmm8, %xmm1			#xmm1 = xmm8
	pxor %xmm9, %xmm1			#xmm1 = xmm8+xmm9
	movdqa 	%xmm12, %xmm7			#xmm7 = xmm12
	pxor %xmm13, %xmm7			#xmm7 = xmm12+xmm13
	pclmulqdq $0x00, %xmm7, %xmm1		#xmm1 = xmm1 * xmm7 = (xmm8+xmm9)*(xmm12+xmm13)
	pxor %xmm0, %xmm1			#xmm1 = xmm1 + xmm0
	pxor %xmm2, %xmm1			#xmm1 = xmm1 + xmm2

	movdqa	%xmm10, %xmm4			#xmm4 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm4		#xmm4 = xmm10 * xmm14
	movdqa	%xmm11, %xmm6			#xmm6 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm6		#xmm6 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm5			#xmm5 = xmm10
	pxor %xmm11, %xmm5			#xmm5 = xmm10 + xmm11
	movdqa 	%xmm14, %xmm7			#xmm7 = xmm14
	pxor %xmm15, %xmm7			#xmm7 = xmm14 + xmm15
	pclmulqdq $0x00, %xmm7, %xmm5		#xmm5 = (xmm10 + xmm11) * (xmm14 + xmm15)
	pxor %xmm4, %xmm5			#xmm5 = xmm5 + xmm4
	pxor %xmm6, %xmm5			#xmm5 = xmm5 + xmm6


	pxor %xmm8, %xmm10			#xmm10 = xmm8 ^ xmm10 = xmm8 + xmm10
	pxor %xmm9, %xmm11			#xmm11 = xmm9 ^ xmm11 = xmm9 + xmm11
	pxor %xmm12, %xmm14			#xmm14 = xmm12 ^ xmm14 = xmm12 + xmm14
	pxor %xmm13, %xmm15			#xmm15 = xmm13 ^ xmm15 = xmm13 + xmm15
	

	movdqa	%xmm10, %xmm7			#xmm7 = xmm10
	pclmulqdq $0x00, %xmm14, %xmm7		#xmm7 = xmm10 * xmm14
	movdqa	%xmm11, %xmm9			#xmm9 = xmm11
	pclmulqdq $0x00, %xmm15, %xmm9		#xmm9 = xmm11 * xmm15
	movdqa 	%xmm10, %xmm8			#xmm8 = xmm10
	pclmulqdq $0x00, %xmm15, %xmm8		#xmm8 = xmm10 * xmm15
	movdqa 	%xmm11, %xmm12			#xmm12 = xmm11
	pclmulqdq $0x00, %xmm14, %xmm12		#xmm12 = xmm11 * xmm14
	pxor %xmm12, %xmm8			#xmm8 = xmm8 ^ xmm12 = xmm10 * xmm15 + xmm11 * xmm14

	pxor %xmm0, %xmm7			#xmm7 = xmm7 ^ xmm0 = xmm7 - xmm0
	pxor %xmm4, %xmm7			#xmm7 = xmm7 ^ xmm4 = xmm7 - xmm4
	pxor %xmm1, %xmm8			#xmm8 = xmm8 ^ xmm1 = xmm8 - xmm1
	pxor %xmm5, %xmm8			#xmm8 = xmm8 ^ xmm5 = xmm8 - xmm5
	pxor %xmm2, %xmm9			#xmm9 = xmm9 ^ xmm2 = xmm9 - xmm2
	pxor %xmm6, %xmm9			#xmm9 = xmm9 ^ xmm6 = xmm9 - xmm6

	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7
	movdqa	%xmm8, %xmm3			#xmm3 = xmm8
	pxor %xmm9, %xmm4			#xmm4 = xmm4 ^ xmm9 = xmm4 + xmm9
	#expand
	movdqa %xmm0, %xmm7			#xmm7 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm1			#xmm1 = xmm1 ^ xmm7 = xmm1 + xmm7

	movdqa %xmm1, %xmm7			#xmm7 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm2			#xmm2 = xmm2 ^ xmm7 = xmm2 + xmm7

	movdqa %xmm2, %xmm7			#xmm7 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7 = xmm3 + xmm7

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask64, %xmm3			#xmm3 = xmm3 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm4			#xmm4 = xmm4 ^ xmm7 = xmm4 + xmm7

	movdqa %xmm4, %xmm7			#xmm7 = xmm4
	pand vMask64, %xmm4			#xmm4 = xmm4 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm5			#xmm5 = xmm5 ^ xmm7 = xmm5 + xmm7

	movdqa %xmm5, %xmm7			#xmm7 = xmm5
	pand vMask64, %xmm5			#xmm5 = xmm5 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64
	pxor %xmm7, %xmm6			#xmm6 = xmm6 ^ xmm7 = xmm6 + xmm7

	movdqa %xmm6, %xmm7			#xmm7 = xmm6
	pand vMask64, %xmm6			#xmm6 = xmm6 & vMask64
	psrldq $8, %xmm7			#xmm7 = xmm7>>64

	#reduce part 1
	pclmulqdq $0x00, residue2, %xmm4	#xmm4 = xmm4 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm5	#xmm5 = xmm5 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm6	#xmm6 = xmm6 * residue2 (=0x12A0)
	pclmulqdq $0x00, residue2, %xmm7	#xmm7 = xmm7 * residue2 (=0x12A0)

	pxor %xmm4, %xmm0			#xmm0 = xmm0 ^ xmm4
	pxor %xmm5, %xmm1			#xmm1 = xmm1 ^ xmm5
	pxor %xmm6, %xmm2			#xmm2 = xmm2 ^ xmm6
	pxor %xmm7, %xmm3			#xmm3 = xmm3 ^ xmm7

	#final reduction
	movdqa %xmm0, %xmm4			#xmm4 = xmm0
	pand vMask64, %xmm0			#xmm0 = xmm0 & vMask64
	psrldq $8, %xmm4			#xmm4 = xmm4 >> 64
	pxor %xmm4, %xmm1			#xmm1 = xmm1 ^ xmm4

	movdqa %xmm1, %xmm5			#xmm5 = xmm1
	pand vMask64, %xmm1			#xmm1 = xmm1 & vMask64
	psrldq $8, %xmm5			#xmm5 = xmm5 >> 64
	pxor %xmm5, %xmm2			#xmm2 = xmm2 ^ xmm5

	movdqa %xmm2, %xmm6			#xmm6 = xmm2
	pand vMask64, %xmm2			#xmm2 = xmm2 & vMask64
	psrldq $8, %xmm6			#xmm6 = xmm6 >> 64
	pxor %xmm6, %xmm3			#xmm3 = xmm3 ^ xmm6

	movdqa %xmm3, %xmm7			#xmm7 = xmm3
	pand vMask59, %xmm3			#xmm3 = xmm3 & vMask59
	psrldq $7, %xmm7			#xmm7 = xmm7 >> 56
	psrlq $3, %xmm7				#xmm7 = xmm7 >> 3 (packed right shift on 64 bit integers by 3)
	pclmulqdq $0x00, residue1, %xmm7	#xmm7 = xmm7 * residue1 (=0x95)
	pxor %xmm7, %xmm0			#xmm0 = xmm0 ^ xmm7

	movdqa	%xmm0,0(%rdx)			
	movdqa	%xmm1,16(%rdx)			
	movdqa	%xmm2,32(%rdx)			
	movdqa	%xmm3,48(%rdx)	

	movq   0(%rsp), %r11
	movq   %r11, %rsp

ret

