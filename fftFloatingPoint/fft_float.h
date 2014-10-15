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

typedef struct SCplxFloatFFT* HCplxFloatFFT;

int CreateCplxFloatFFT(HCplxFloatFFT* h, int numPoints);

int ComputeCplxFloatFFT(HCplxFloatFFT h, const float* in, float* out, EDirection direction);

void DisposeCplxFloatFFT(HCplxFloatFFT h);


typedef HCplxFloatFFT HRealFloatFFT;

int CreateRealFloatFFT(HRealFloatFFT* h, int numPoints);

int ComputeRealFloatFFT(HRealFloatFFT h, const float* in, float* out, EDirection direction);

void DisposeRealFloatFFT(HRealFloatFFT h);

void ComputeFloatEnergies(const float* spec, float* energies, int size);

#ifdef __cplusplus
}
#endif
  
#endif
