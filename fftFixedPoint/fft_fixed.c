/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include <stdlib.h>
#include "functions.h"
#include "fft_fixed.h"
#include "fft_fixed_impl.h"

struct SCplxFixedFFT
{
	int numPoints;
	int scale;
	PackedType* cosSinTable;
	int* tempBuffer;
} SCplxFixedFFT;


int CreateCplxFixedFFT(HCplxFixedFFT* h, int numPoints)
{
	int numBits = 16;
	*h = (HCplxFixedFFT)malloc(sizeof(SCplxFixedFFT));
	if(NULL == *h)
		return -1;
	(*h)->tempBuffer = NULL;

	(*h)->numPoints   = numPoints;
	(*h)->cosSinTable = CreateFixedCosSinTable(numPoints, numBits);
	(*h)->tempBuffer = (int*) malloc (numPoints *sizeof(int));
	if ((*h)->tempBuffer == NULL)
		return -1;
	(*h)->scale = GetFixedScalingBits(numPoints);
	
	return 0;
}


int ComputeCplxFixedFFT(HCplxFixedFFT h, 
                   const int* in, 
                   int* out, 
                   EDirection dir)
{
	Radix4_CplxFixedFFT(in, out, h->cosSinTable, 1, h->numPoints);
	
	if (dir == Inverse)
	{
		InvertFixedOutput(out, h->numPoints);
	}
	
	return 0;
}


void DisposeCplxFixedFFT(HCplxFixedFFT h)
{
	if(h)
	{
		DisposeFixedCosSinTable(h->cosSinTable);
		if (h->tempBuffer)
			free(h->tempBuffer);
		h->tempBuffer = NULL;
		free(h);
	}
}



int CreateRealFixedFFT(HRealFixedFFT* h, int numPoints)
{
	return CreateCplxFixedFFT(h, numPoints);
}


int ComputeRealFixedFFT(HRealFixedFFT h, 
                    const int* in, 
                    int* out, 
                    EDirection dir)
{
	int i;
	if (dir == Forward)
	{
		Radix4_CplxFixedFFT(in, out, h->cosSinTable, 2, h->numPoints/2);
		CplxFixedToReal(out, h->cosSinTable, h->numPoints/2);
	}
	else
	{
		ComplexFixedToRealInversion(in, h->tempBuffer, h->cosSinTable, h->numPoints/2);
		// Complex conjugate of the input: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			h->tempBuffer[i] = -h->tempBuffer[i];
		Radix4_CplxFixedFFT(h->tempBuffer, out, h->cosSinTable, 2, h->numPoints/2);
		// Complex conjugate of the output: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			out[i] = -out[i];
		// Scale the output - the 1/FFTSize scaling occurs in the ComplexToRealInversion for floating point
		for (i=0; i<h->numPoints; i++)
			out[i] = out[i]<<h->scale;		
	}
	
	return 0;
}

void DisposeRealFixedFFT(HRealFixedFFT h)
{
	DisposeCplxFixedFFT(h);
}

void ComputeFixedEnergies(const int* spec, int* energies, int N)
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
