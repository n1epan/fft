#ifndef __FFT_BFLY_FIX_H__
#define __FFT_BFLY_FIX_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include "functions.h"

#define ScaleFxpRight(x,s) (x>>s)
#define ScaleFxpLeft(x,s)  (x<<s)

/* Basic Radix 2 W0 butterfly macro */
#define Radix2_Bfly_W0(x0_r, x0_i, x1_r, x1_i, s) \
{	\
	Type t; \
	x0_r >>= 1; \
	t = x0_r + (x1_r >> s); \
	x1_r = x0_r - (x1_r >> s); \
	x0_r = t; \
	x0_i >>= 1; \
	t = x0_i + (x1_i >> s); \
	x1_i = x0_i - (x1_i >> s); \
	x0_i = t; \
}

/* Basic Radix 4 W0 butterfly macro */
#define Radix4_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, s) \
{	\
  Type t; \
  x2_r >>= s; \
  x2_i >>= s; \
  x2_r = x2_r + (x3_r >> s); \
  x2_i = x2_i + (x3_i >> s); \
  x3_r = x2_r - (x3_r >> (s-1)); \
  x3_i = x2_i - (x3_i >> (s-1)); \
\
  x0_r >>= 2; \
  x0_i >>= 2; \
  x0_r = x0_r + (x1_r >> s); \
  x0_i = x0_i + (x1_i >> s); \
  x1_r = x0_r - (x1_r >> (s-1)); \
  x1_i = x0_i - (x1_i >> (s-1)); \
\
  x0_r = x0_r + x2_r; \
  x0_i = x0_i + x2_i; \
  x2_r = x0_r - (x2_r << 1); \
  x2_i = x0_i - (x2_i << 1); \
\
  t = x3_r; \
  x1_r = x1_r + x3_i; \
  x1_i = x1_i - x3_r; \
  x3_r = x1_r - (x3_i << 1); \
  x3_i = x1_i + (t << 1); \
} \


/* 
Complex mutliplication macro:
 d = a*b with b packed type
 d_r = a_r*b_r - a_i*b_i
 d_i = a_i*b_r + a_r*b_i
 */
#define CplxMul(d_r, d_i, a_r, a_i, b) \
{ \
  d_r = SmulWLo_SW_SL(a_i, b); \
  d_i = SmulWLo_SW_SL(a_r, b); \
  d_r = 0 - d_r; \
  d_r = SmulAddWHi_SW_SL(a_r, b, d_r); \
  d_i = SmulAddWHi_SW_SL(a_i, b, d_i); \
}


/* 
Complex conjugate mutliplication macro:
 d = a*(conj)b with b packed type
 d_r = a_r*b_r + a_i*b_i
 d_i = a_i*b_r - a_r*b_i
 */
#define CplxMulConj(d_r, d_i, a_r, a_i, b) \
{ \
  d_r = SmulWHi_SW_SL(a_r, b); \
  d_i = SmulWLo_SW_SL(a_r, b); \
  d_r = SmulAddWLo_SW_SL(a_i, b, d_r); \
  d_i = 0 - d_i; \
  d_i = SmulAddWHi_SW_SL(a_i, b, d_i); \
}

#endif
