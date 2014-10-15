/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include <stdlib.h>
#include "functions.h"
#include "fft.h"
#include "fft_impl.h"

struct SCplxFFT
{
	int numPoints;
	int scale;
	PackedType* cosSinTable;
	Type* tempBuffer;
} SCplxFFT;


int CreateCplxFFT(HCplxFFT* h, int numPoints)
{
	int numBits = 16;
	*h = (HCplxFFT)malloc(sizeof(SCplxFFT));
	if(NULL == *h)
		return -1;
	(*h)->tempBuffer = NULL;

	(*h)->numPoints   = numPoints;
	(*h)->cosSinTable = CreateCosSinTable(numPoints, numBits);
	(*h)->tempBuffer = (Type*) malloc (numPoints *sizeof(Type));
	if ((*h)->tempBuffer == NULL)
		return -1;
	(*h)->scale = GetScalingBits(numPoints);
	
	return 0;
}


int ComputeCplxFFT(HCplxFFT h, 
                   const Type* in, 
                   Type* out, 
                   EDirection dir)
{
	Radix4_CplxFFT(in, out, h->cosSinTable, 1, h->numPoints);
	
	if (dir == Inverse)
	{
		InvertOutput(out, h->numPoints);
	}
	
	return 0;
}


void DisposeCplxFFT(HCplxFFT h)
{
	if(h)
	{
		DisposeCosSinTable(h->cosSinTable);
		if (h->tempBuffer)
			free(h->tempBuffer);
		h->tempBuffer = NULL;
		free(h);
	}
}



int CreateRealFFT(HRealFFT* h, int numPoints)
{
	return CreateCplxFFT(h, numPoints);
}


int ComputeRealFFT(HRealFFT h, 
                    const Type* in, 
                    Type* out, 
                    EDirection dir)
{
	int i;
	if (dir == Forward)
	{
		Radix4_CplxFFT(in, out, h->cosSinTable, 2, h->numPoints/2);
		CplxToReal(out, h->cosSinTable, h->numPoints/2);
	}
	else
	{
		ComplexToRealInversion(in, h->tempBuffer, h->cosSinTable, h->numPoints/2, h->scale);
		// Complex conjugate of the input: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			h->tempBuffer[i] = -h->tempBuffer[i];
		Radix4_CplxFFT(h->tempBuffer, out, h->cosSinTable, 2, h->numPoints/2);
		// The RotateSpectrum function could be used instead of the two complex conjugation loops
//		RotateSpectrum(in,h->numPoints/2);
		// Complex conjugate of the output: change the sign on the odd samples(imaginary samples)
		for (i=1; i<h->numPoints; i=i+2)
			out[i] = -out[i];
#ifdef FXP_32_BIT
		// Scale the output - the 1/FFTSize scaling occurs in the ComplexToRealInversion for floating point
		for (i=0; i<h->numPoints; i++)
			out[i] = out[i]<<h->scale;
#endif
		

	}
	
	return 0;
}

void DisposeRealFFT(HRealFFT h)
{
	DisposeCplxFFT(h);
}

void ComputeEnergies(const Type* spec, Type* energies, int N)
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
