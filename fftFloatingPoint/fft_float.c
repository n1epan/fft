/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include <stdlib.h>
#include "functions.h"
#include "fft_float.h"
#include "fft_float_impl.h"

struct SCplxFloatFFT
{
	int numPoints;
	int scale;
	PackedType* cosSinTable;
	float* tempBuffer;
} SCplxFloatFFT;


int CreateCplxFloatFFT(HCplxFloatFFT* h, int numPoints)
{
	int numBits = 16;
	*h = (HCplxFloatFFT)malloc(sizeof(SCplxFloatFFT));
	if(NULL == *h)
		return -1;
	(*h)->tempBuffer = NULL;

	(*h)->numPoints   = numPoints;
	(*h)->cosSinTable = CreateFloatCosSinTable(numPoints, numBits);
	(*h)->tempBuffer = (float*) malloc (numPoints *sizeof(float));
	if ((*h)->tempBuffer == NULL)
		return -1;
	(*h)->scale = GetFloatScalingBits(numPoints);
	
	return 0;
}


int ComputeCplxFloatFFT(HCplxFloatFFT h, 
                   const float* in, 
                   float* out, 
                   EDirection dir)
{
	Radix4_CplxFloatFFT(in, out, h->cosSinTable, 1, h->numPoints);
	
	if (dir == Inverse)
	{
		InvertFloatOutput(out, h->numPoints);
	}
	
	return 0;
}


void DisposeCplxFloatFFT(HCplxFloatFFT h)
{
	if(h)
	{
		DisposeFloatCosSinTable(h->cosSinTable);
		if (h->tempBuffer)
			free(h->tempBuffer);
		h->tempBuffer = NULL;
		free(h);
	}
}



int CreateRealFloatFFT(HRealFloatFFT* h, int numPoints)
{
	return CreateCplxFloatFFT(h, numPoints);
}


int ComputeRealFloatFFT(HRealFloatFFT h, 
                    const float* in, 
                    float* out, 
                    EDirection dir)
{
	int i;
	if (dir == Forward)
	{
		Radix4_CplxFloatFFT(in, out, h->cosSinTable, 2, h->numPoints/2);
		CplxToFloatReal(out, h->cosSinTable, h->numPoints/2);
	}
	else
	{
		ComplexToFloatRealInversion(in, h->tempBuffer, h->cosSinTable, h->numPoints/2);
		// Complex conjugate of the input: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			h->tempBuffer[i] = -h->tempBuffer[i];
		Radix4_CplxFloatFFT(h->tempBuffer, out, h->cosSinTable, 2, h->numPoints/2);
		// Complex conjugate of the output: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			out[i] = -out[i];		
	}
	
	return 0;
}

void DisposeRealFloatFFT(HRealFloatFFT h)
{
	DisposeCplxFloatFFT(h);
}

void ComputeFloatEnergies(const float* spec, float* energies, int N)
{
	int specIndex = 0;
	int enIndex   = 1;
	
	energies[0]   = spec[0]*spec[0];
	energies[N/2] = spec[1]*spec[1];  //energy 513 discarded?

	for(specIndex = 2; specIndex < N ; specIndex+=2 )
	{
		energies[enIndex++] = spec[specIndex]*spec[specIndex] + spec[specIndex + 1]*spec[specIndex + 1];
		
	}
	
}
