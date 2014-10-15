#ifndef __FFT_IMPL_H__
#define __FFT_IMPL_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include "fft.h"


#ifdef FXP_32_BIT                                                              
typedef Type PackedType;
#else
typedef struct PackedType 
{
	float x; /* Top value */
	float y; /* Bottom value */
} PackedType;
#endif


extern PackedType* CreateCosSinTable(int N_FFT, int numBits);

extern int GetScalingBits(int N_FFT);

extern void DisposeCosSinTable(PackedType* cosSinTable);

extern int Radix4_CplxFFT(const Type* x, Type* y,
                          const PackedType* cosSinTable,
                          int strides, int N_FFT);

extern void CplxToReal(Type* x, const PackedType* cosSinTable, int N_FFT);

extern void ComplexToRealInversion(const Type* in, Type* out, const PackedType* cosSinTable, int N_FFT);

/*extern void RealInverseFFT(const Type* input, Type* out, const PackedType* cosSinTable, int N_FFT);*/

/*extern void RotateSpectrum (Type* x, int N_FFT);*/

#endif
