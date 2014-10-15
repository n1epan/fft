#ifndef __CPLXFFT_H__
#define __CPLXFFT_H__

/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#ifdef __cplusplus
extern "C" {
#endif

#define Type int

enum EDirection
{
	Forward = 0,
	Inverse = 1
};

typedef enum EDirection EDirection;

typedef struct SCplxFixedFFT* HCplxFixedFFT;

int CreateCplxFixedFFT(HCplxFixedFFT* h, int numPoints);

int ComputeCplxFixedFFT(HCplxFixedFFT h, const int* in, int* out, EDirection direction);

void DisposeCplxFixedFFT(HCplxFixedFFT h);


typedef HCplxFixedFFT HRealFixedFFT;

int CreateRealFixedFFT(HRealFixedFFT* h, int numPoints);

int ComputeRealFixedFFT(HRealFixedFFT h, const int* in, int* out, EDirection direction);

void DisposeRealFixedFFT(HRealFixedFFT h);

void ComputeFixedEnergies(const int* spec, int* energies, int size);

#ifdef __cplusplus
}
#endif
  
#endif
