#ifndef __FFT_BFLY_FLOAT_H__
#define __FFT_BFLY_FLOAT_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

/*
 * Radix 4 butterfly for float implementation
 */

/* Basic Radix 2 W0 butterfly macro */
#define Radix2_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, s) \
{	\
	float t; \
	t    = x0_r + x1_r; \
	x1_r = x0_r - x1_r; \
	x0_r = t; \
	t    = x0_i + x1_i; \
	x1_i = x0_i - x1_i; \
	x0_i = t; \
}


/* Basic Radix 4 W0 butterfly macro */
#define Radix4_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, s) \
{	\
	float t; \
	t    = x2_r + x3_r; \
	x3_r = x2_r - x3_r; \
	x2_r = t; \
	t    = x2_i + x3_i; \
	x3_i = x2_i - x3_i; \
	x2_i = t; \
\
	t    = x0_r + x1_r; \
  x1_r = x0_r - x1_r; \
  x0_r = t; \
  t    = x0_i + x1_i; \
  x1_i = x0_i - x1_i; \
  x0_i = t; \
\
  t    = x0_r + x2_r; \
  x2_r = x0_r - x2_r; \
  x0_r = t; \
  t    = x0_i + x2_i; \
  x2_i = x0_i - x2_i; \
  x0_i = t; \
\
	t    = x3_i; \
	x3_i = x1_i + x3_r; \
	x1_i = x1_i - x3_r; \
	x3_r = x1_r - t; \
	x1_r = x1_r + t; \
} \


/* 
Complex mutliplication macro:
 d = a*b with b packed type
 d_r = a_r*b_r - a_i*b_i
 d_i = a_i*b_r + a_r*b_i
 */
#define CplxMulFloat(d_r, d_i, a_r, a_i, b) \
{ \
	d_r = a_r * b.x - a_i * b.y; \
	d_i = a_i * b.x + a_r * b.y; \
}


/* 
Complex conjugate mutliplication macro:
 d = a*(conj)b with b packed type
 d_r = a_r*b_r + a_i*b_i
 d_i = a_i*b_r - a_r*b_i
 */
#define CplxMulConjFloat(d_r, d_i, a_r, a_i, b) \
{ \
	d_r = a_r * b.x + a_i * b.y; \
	d_i = a_i * b.x - a_r * b.y; \
}


#endif
