#ifndef __FFT_IMPL_H__
#define __FFT_IMPL_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include "fft_float.h"


typedef struct PackedType 
{
	float x; /* Top value */
	float y; /* Bottom value */
} PackedType;

extern PackedType* CreateFloatCosSinTable(int N_FFT, int numBits);

extern int GetFloatScalingBits(int N_FFT);

extern void DisposeFloatCosSinTable(PackedType* cosSinTable);

extern int Radix4_CplxFloatFFT(const float* x, float* y,
                          const PackedType* cosSinTable,
                          int strides, int N_FFT);

extern void CplxToFloatReal(float* x, const PackedType* cosSinTable, int N_FFT);

extern void ComplexToFloatRealInversion(const float* in, float* out, const PackedType* cosSinTable, int N_FFT);

extern void InvertFloatOutput(float* x, int N_FFT);


#endif
