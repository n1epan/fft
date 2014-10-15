#ifndef __FFT_IMPL_H__
#define __FFT_IMPL_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include "fft_fixed.h"

                                                         
typedef int PackedType;


extern PackedType* CreateFixedCosSinTable(int N_FFT, int numBits);

extern int GetFixedScalingBits(int N_FFT);

extern void DisposeFixedCosSinTable(PackedType* cosSinTable);

extern int Radix4_CplxFixedFFT(const int* x, int* y,
                          const PackedType* cosSinTable,
                          int strides, int N_FFT);

extern void CplxFixedToReal(int* x, const PackedType* cosSinTable, int N_FFT);

extern void ComplexFixedToRealInversion(const int* in, int* out, const PackedType* cosSinTable, int N_FFT);

extern void InvertFixedOutput(int* x, int N_FFT);


#endif
