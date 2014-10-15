/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft_float_impl.h"
#include "functions.h"
#include "fft_bfly_float.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

FILE* g_out = NULL;


PackedType* CreateFloatCosSinTable(int numPoints, int numBits)
{
  int i;
  
  int size = 3*(numPoints/4 - 1) + 1;
  
  PackedType* cosSinTable = (PackedType*)malloc(sizeof(PackedType)*size);
  
  if(NULL == cosSinTable)
    return NULL;
  
  for(i = 0; i < size; ++i)
  {
    cosSinTable[i].x = cos(2*M_PI*(float)i/numPoints);
    cosSinTable[i].y = sin(2*M_PI*(float)i/numPoints);
  }

  return cosSinTable;
}


void DisposeFloatCosSinTable(PackedType* cosSinTable)
{
  if(cosSinTable)
    free(cosSinTable);
}

int GetFloatScalingBits(int N_FFT)
{
	int count = 0;
	while (N_FFT>2)
	{
		N_FFT = N_FFT>>1;
		count++;
	}
	return count;
}


int Radix4_CplxFloatFFT(const float* x, float* y,
                   const PackedType* cosSinTable,
                   int trigTableStrides, int N_FFT)
{
  int idx;
  int span;
  int twiddle;
  int strides;
	
  float x0_r, x0_i;
  float x1_r, x1_i;
  float x2_r, x2_i;
  float x3_r, x3_i;
  
  
  /*
   First stage:
   load data in reversed carry arithmetic (modulo 4) 
   and compute multiply free butterflies. Results are stored
   in natural order.
   */
  
  idx = span = 0;

  do
  {
    /* Load */
    x0_r = x[2*span];
    x0_i = x[2*span+1];
    span += N_FFT/4;
    x2_r = x[2*span];
    x2_i = x[2*span+1];
    span += N_FFT/4;    
    x1_r = x[2*span];
    x1_i = x[2*span+1];
    span += N_FFT/4;
    x3_r = x[2*span];
    x3_i = x[2*span+1];

    /* W0 butterfly */
    Radix4_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);
    
    /* Store */
    y[2*idx]   = x0_r;
    y[2*idx+1] = x0_i;
    idx += 1;
    y[2*idx]   = x1_r;
    y[2*idx+1] = x1_i;
    idx += 1;
    y[2*idx]   = x2_r;
    y[2*idx+1] = x2_i;
    idx += 1;
    y[2*idx]   = x3_r;
    y[2*idx+1] = x3_i;
    idx += 1;
   
    /*
     Reverse carry modulo N/4:
     should be span = N_FFT/4 - span -1
     but since span was increased by 3*N_FFT/4,
     the above formula reduces to span = N_FFT - span -1
    */
    span = N_FFT - span - 1;
    span = (span ^ (0x7fffffffUL >> CountLeadingZeros(span)));

  } while (idx < N_FFT);
  
  span    = 4;
  strides = trigTableStrides*N_FFT/16;

  while(span < N_FFT/2)
  {
    /* Multiplication free butterflies */
    idx = 3*span;

    do
    {
      /* Load */
      x3_r = y[2*idx];
      x3_i = y[2*idx+1];
      idx -= span;
      x2_r = y[2*idx];
      x2_i = y[2*idx+1];
      idx -= span;
      x1_r = y[2*idx];
      x1_i = y[2*idx+1];
      idx -= span;
      x0_r = y[2*idx];
      x0_i = y[2*idx+1];

      /* W0 butterfly */
      Radix4_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);

      /* Store */
      y[2*idx]   = x0_r;
      y[2*idx+1] = x0_i;
      idx += span;
      y[2*idx]   = x1_r;
      y[2*idx+1] = x1_i;
      idx += span;
      y[2*idx]   = x2_r;
      y[2*idx+1] = x2_i;
      idx += span;
      y[2*idx]   = x3_r;
      y[2*idx+1] = x3_i;
      idx += 4*span;
      
    } while (idx < N_FFT);
    
    twiddle = 1;

    do
    {
      PackedType cosSin1, cosSin2, cosSin3;
      
      /* Load packed cos, sin coefficients */
      cosSin1 = cosSinTable[strides*twiddle];
      cosSin2 = cosSinTable[2*strides*twiddle];
      cosSin3 = cosSinTable[3*strides*twiddle];

      /* 
        Advance data index for subsequent 
        reverse reading
       */
      idx = 3*span+twiddle;

      do
      {
        x2_r = y[2*idx];
        x2_i = y[2*idx+1];
        idx -= span;
        x1_r = y[2*idx];
        x1_i = y[2*idx+1];
        idx -= span;
        x0_r = y[2*idx];
        x0_i = y[2*idx+1];
        idx -= span;
        
        /* Muliplication with twiddle factors */
        CplxMulConjFloat(x3_r, x3_i, x2_r, x2_i, cosSin3);
        CplxMulConjFloat(x2_r, x2_i, x1_r, x1_i, cosSin1);
        CplxMulConjFloat(x1_r, x1_i, x0_r, x0_i, cosSin2);
        
        x0_r = y[2*idx];
        x0_i = y[2*idx+1];
            
        /* W0 butterfly */
        Radix4_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 1);
        
        y[2*idx]   = x0_r;
        y[2*idx+1] = x0_i;
        idx += span;
        y[2*idx]   = x1_r;
        y[2*idx+1] = x1_i;
        idx += span;
        y[2*idx]   = x2_r;
        y[2*idx+1] = x2_i;
        idx += span;
        y[2*idx]   = x3_r;
        y[2*idx+1] = x3_i;
        idx += 4*span;
      
      } while (idx < N_FFT);
      
      ++twiddle;
      
    } while(twiddle < span);
    strides /= 4;
    span    *= 4;
  }
  
	
	/* 
		If N_FFT is exressed as 4^N * 2, perform a last radix2 stage
	 */
	if(span < N_FFT)
	{
		idx = span;
		
		x1_r = y[2*idx];
		x1_i = y[2*idx+1];
		idx -= span;
		x0_r = y[2*idx];
		x0_i = y[2*idx+1];
				
		/* W0 radix 2 butterfly */
		Radix2_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 1);
		
		y[2*idx]   = x0_r;
		y[2*idx+1] = x0_i;
		idx += span;
		y[2*idx]   = x1_r;
		y[2*idx+1] = x1_i;
		
		twiddle = 1;
		strides = trigTableStrides;
		
		do
		{
			PackedType cosSin1 = cosSinTable[strides*twiddle];
			
			++idx;
			
			x0_r = y[2*idx];
			x0_i = y[2*idx+1];
			idx -= span;
			
			CplxMulConjFloat(x1_r, x1_i, x0_r, x0_i, cosSin1);

			x0_r = y[2*idx];
			x0_i = y[2*idx+1];
						
			/* W0 radix 2 butterfly */
			Radix2_Float_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 0);
			
			y[2*idx]   = x0_r;
			y[2*idx+1] = x0_i;
			idx += span;
			y[2*idx]   = x1_r;
			y[2*idx+1] = x1_i;
			
			++twiddle;
			
		} while(twiddle < span);
	}
	
  return 0;
}

// Since the even samples go in as the real input and the odd samples go in as the imaginary input
// to a complex FFT, this function is needed to remap the output, Z(k), to the correct F(k).
// So, Feven(k) = (Z(k) + conj Z(N/2 - k))*0.5
// Fodd(k) = -j(Z(k) - conj Z(N/2 - k))*0.5
// Now combine these two using a Decimation In Time butterfly:
// F(k) = Feven(k) + e^-j2pik/N*Fodd(k)
// F(k) = 1/2*(Z(k) + conj Z(N/2 -k) - je^-j2pik/N(Z(k) - conj Z(N/2 - k))
void CplxToFloatReal(float* x, const PackedType* cosSinTable, int N_FFT)
{
  int idx;
  
  float x0_r, x0_i;
  float x1_r, x1_i;
  float x2_r, x2_i;
  PackedType cosSin;
  
  /*
   X(0) = Re[X(0)] + Im[X(0)]
   X(N) = Re[X(0)] - Im[X(0)]
   Since X(0) and X(N) are real, we map X(N) to
   Im[X(N)]
   */
  x0_r = x[0];
  x0_i = x[1];
  
  x1_r = x0_r + x0_i;
  x0_i = x0_r - x0_i;
	x0_r = x1_r;
  
  x[0] = x0_r;
  x[1] = x0_i;
  
  
  idx = 1;
  
  do
  {
    x0_r = x[2*idx];
    x0_i = x[2*idx+1];
    x1_r = x[2*(N_FFT-idx)];
    x1_i = x[2*(N_FFT-idx)+1];

    /*
     x0 = 0.5*(X(n) + X'(N-n))
     x1 = 0.5*(X(n) - X'(N-n))
     */
    x0_r = x0_r/(Type)2;
		x0_i = x0_i/(Type)2;
		x0_r = x0_r + x1_r/(Type)2;
    x0_i = x0_i - x1_i/(Type)2;
    x1_r = x0_r - x1_r;
    x1_i = x0_i + x1_i;
    
    cosSin = cosSinTable[idx];
    
    CplxMulConjFloat(x2_r, x2_i, x1_r, x1_i, cosSin);
    
    x1_r = x0_r + x2_i;
    x1_i = x0_i - x2_r;
    x0_r = x0_r - x2_i;
    x0_i = -(x0_i + x2_r);
    
    x[2*idx]           = x1_r;
    x[2*idx+1]         = x1_i;
    x[2*(N_FFT-idx)]   = x0_r;
    x[2*(N_FFT-idx)+1] = x0_i;
    ++idx;
    
  } while(idx < N_FFT/2);
  
  /* Value at N_FFT/2: only change sign of imaginary part */
  x[2*idx]   = x[2*idx];
  x[2*idx+1] = -x[2*idx+1];

}

// This function performs the inverse of the above complex to real function
// This gets us back to a single complex spectrum to act as an input to the complex FFT
// The output of the real FFT is F(k), which is made up of the Feven(k) and Fodd(k) values
// To get these values back perform the inverse matrix multiplication of the above function
// which yields:
// Feven(k) = 1/2*(F(k) + conjF(N/2-k))
// Fodd(k) = 1/2*(F(k) - conjF(N/2-k))*e^j2pik/N
// This is then combined to give us the Complex(k) that the complex FFT originally output.
// Complex(k) = Feven(k) + j*Fodd(k)
// This function also scales the samples to 1/N since this is used for the inverse FFT
// In place implementation
void ComplexToFloatRealInversion(const float* in, float* out, const PackedType* cosSinTable, int N_FFT)
{
	int idx, N_HALF;
	float x0_r, x0_i, x1_r, x1_i, x2_r, x2_i;
	float diag,onePlusSin,negOneMinSin,cosine;
	PackedType cosSin, cosVal, sinVal;
	float scale;

	N_HALF = N_FFT>>1;
	
	scale = 1.0f/(float)2.0f*N_FFT;
	out[0] = (in[0] + in[1])*scale*0.5f;
	out[1] = (in[0] - in[1])*scale*0.5f;
	out[N_FFT] = in[N_FFT]*scale;
	out[N_FFT+1] = -in[N_FFT+1]*scale;

	idx = 1;
  
  do
  {
	x0_r = in[2*idx];
	x0_i = in[2*idx+1];
	x1_r = in[2*(N_FFT-idx)];
	x1_i = in[2*(N_FFT-idx)+1];

	// floating point - inverse matrix for the CplxToReal function
    cosSin = cosSinTable[idx];
	diag = (1.0f - cosSin.y)*0.5f;
	onePlusSin = (1.0f + cosSin.y)*0.5f;
	negOneMinSin = (-1.0f - cosSin.y)*0.5f;
	cosine = cosSin.x*0.5f;

    out[2*idx]           = diag*x0_r - cosine*x0_i + onePlusSin*x1_r - cosine*x1_i;
    out[2*idx+1]         = cosine*x0_r + diag*x0_i - cosine*x1_r + negOneMinSin*x1_i;
    out[2*(N_FFT-idx)]   = onePlusSin*x0_r + cosine*x0_i + diag*x1_r + cosine*x1_i;
    out[2*(N_FFT-idx)+1] = cosine*x0_r + negOneMinSin*x0_i - cosine*x1_r + diag*x1_i;
	
	out[2*idx]			*= scale;
	out[2*idx+1]		*= scale;
	out[2*(N_FFT-idx)]	*= scale;
	out[2*(N_FFT-idx)+1]*= scale; 

    ++idx;
  } while(idx < N_HALF);
}

void InvertFloatOutput(float* x, int N_FFT)
{
	int i;
	float scale = 1.0f / N_FFT;
	
	/* the first pair are merely scaled */
	x[0] = x[0] * scale;
	x[1] = x[1] * scale;
	
	/* the rest are scaled and swapped about the middle */
	for (i = 1; i <= N_FFT / 2; i++)
	{
		int a_index = 2 * i;
		int b_index = 2 * (N_FFT - i);
		
		float a_real = x[a_index];
		float a_imag = x[a_index + 1];
		float b_real = x[b_index];
		float b_imag = x[b_index + 1];
		
		x[a_index]     = (float)(b_real * scale);
		x[a_index + 1] = (float)(b_imag * scale);
		x[b_index]     = (float)(a_real * scale);
		x[b_index + 1] = (float)(a_imag * scale);
	}
}



