#ifndef __CPLXFFT_H__
#define __CPLXFFT_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#ifdef FXP_32_BIT

/* 32 bit integer data type */
#define Type int
#else

/* Float data type */
#define Type float
#endif


#ifdef __cplusplus
extern "C" {
#endif

enum EDirection
{
	Forward = 0,
	Inverse = 1
};

typedef enum EDirection EDirection;

typedef struct SCplxFFT* HCplxFFT;

int CreateCplxFFT(HCplxFFT* h, int numPoints);

int ComputeCplxFFT(HCplxFFT h, const Type* in, Type* out, EDirection direction);

void DisposeCplxFFT(HCplxFFT h);


typedef HCplxFFT HRealFFT;

int CreateRealFFT(HRealFFT* h, int numPoints);

int ComputeRealFFT(HRealFFT h, const Type* in, Type* out, EDirection direction);

void DisposeRealFFT(HRealFFT h);

void ComputeEnergies(const Type* spec, Type* energies, int size);

#ifdef __cplusplus
}
#endif
  
#endif
