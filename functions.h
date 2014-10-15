#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#if defined(_WIN32_WCE) && defined(ARM)
#include <cmnintrin.h>
#include <armintr.h>

#define USE_INLINE

#define CountLeadingZeros(a) _CountLeadingZeros(a)
#define SmulWLo_SW_SL(a,b) _SmulWLo_SW_SL(a,b)
#define SmulWHi_SW_SL(a,b) _SmulWHi_SW_SL(a,b)

/* 
  Note: the additive term "c" must be past as the first argument
  to the M$ ARM-intrinsics implementation, whereas the ARM assembler
  instruction takes it at the last argument
 */
#define SmulAddWLo_SW_SL(a,b,c) _SmulAddWLo_SW_SL(c,a,b)
#define SmulAddWHi_SW_SL(a,b,c) _SmulAddWHi_SW_SL(c,a,b)

#elif defined(__GNUC__) && defined(__arm__)
#define USE_INLINE

inline int SmulWLo_SW_SL(int a, int b)
{
	int r;
	asm volatile("smulwb %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
	return r;
}

inline int SmulWHi_SW_SL(int a, int b)
{
	int r;
	asm volatile("smulwt %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
	return r;
}

inline int SmulAddWLo_SW_SL(int a, int b, int c)
{
	int r;
	asm volatile("smlawb %0, %1, %2, %3" : "=r"(r) : "r"(a), "r"(b), "r"(c));
	return r;
}

inline int SmulAddWHi_SW_SL(int a, int b, int c)
{
	int r;
	asm volatile("smlawt %0, %1, %2, %3" : "=r"(r) : "r"(a), "r"(b), "r"(c));
	return r;
}

inline int CountLeadingZeros(int a)
{
	int r;
	asm volatile("clz %0, %1" : "=r"(r) : "r"(a));
	return r;
}


/* 
  Generic C function declarations
 */
#else  /*  #elif defined(__GNUC__) && defined(__arm__) */

extern int CountLeadingZeros(int x);

extern int FixedPointMul(int a, int b);

extern int FixedPointMulAndAdd(int a, int b, int c);

extern int FixedPointNegMulAndAdd(int a, int b, int c);

extern int SmulWLo_SW_SL(int a, int b);

extern int SmulWHi_SW_SL(int a, int b);

extern int SmulAddWLo_SW_SL(int a, int b, int c);

extern int SmulAddWHi_SW_SL(int a, int b, int c);

#endif

extern int RoundFloatToInt(float x, int numBits);

#endif
