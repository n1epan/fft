/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include "functions.h"


/*
 Functions which are already declared as inline/ intrinsics
 */
#ifndef USE_INLINE
int CountLeadingZeros(int x)
{
	int count = 32;
	while(x)
	{
		x = ((unsigned int)x >> 1);
		--count;
	}
	return count;
}

int FixedPointMul(int a, int b)
{
	int r;
	r = (a >> 16) * (b >> 16); // hi * hi goes to hi bits
	r += ((a & 0x0000ffff) * (b >> 16)) >> 16; // hi * lo goes to lo bits
	r += ((a >> 16) * (b & 0x0000ffff)) >> 16; // hi * lo add to lo bits
	// skip the lo * lo since it goes into the bits below zero :)
	return r;
}

int FixedPointMulAndAdd(int a, int b, int c)
{
	// this does not check for over flow
	return (FixedPointMul(a, b) + c);
}

int FixedPointNegMulAndAdd(int a, int b, int c)
{
	// this does not check for over flow
	// subtract the a*b product from c
	return (c - FixedPointMul(a,b));
}

int SmulWLo_SW_SL(int a, int b)
{
	int r;
	int bs = (b << 16) >> 16;
	r  = (a >> 16) * bs;
	r += ((a & 0X0000FFFF) * bs) >> 16;
	
	return r;
}


int SmulWHi_SW_SL(int a, int b)
{
	int r;
	int bs = b >> 16;
	r  = (a >> 16) * bs;
	r += ((a & 0X0000FFFF) * bs) >> 16;
	
	return r;
}


int SmulAddWLo_SW_SL(int a, int b, int c)
{
	return (SmulWLo_SW_SL(a, b) + c);
}


int SmulAddWHi_SW_SL(int a, int b, int c)
{
	return (SmulWHi_SW_SL(a, b) + c);
}

#endif  /* #ifndef USE_INLINE */

/*
 Common, non inlined functions
 */
int RoundFloatToInt (float x, int numBits)
{
	int y;

	/* Positive float values */
	if (x >= 0.f)
	{
		x += 0.5f;
		if (x >= (unsigned int)(1<<numBits))
			y = (1<<numBits)-1;
		else
			y = (int)x;    
	}
	/* Positive float values */ 
	else
	{
		x -= 0.5f;
		if ( x< ((-1)<<numBits))
			y = ((-1)<<numBits)+1;
		else
			y = (int)x;    
	}
	return y;
}
