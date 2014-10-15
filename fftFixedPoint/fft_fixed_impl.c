/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */


#include <stdlib.h>
#include <math.h>
#include "fft_fixed_impl.h"
#include "functions.h"
#include "fft_bfly_fix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif



PackedType* CreateFixedCosSinTable(int numPoints, int numBits)
{
  const int mask = (1<<(numBits))-1;
  float scale = (float)pow(2.0, numBits-1);
  int i;
  
  int size = 3*(numPoints/4 - 1) + 1;
  
  PackedType* cosSinTable = (PackedType*)malloc(sizeof(PackedType)*size);
  
  if(NULL == cosSinTable)
    return NULL;
  
  for(i = 0; i < size; ++i)
  {
    cosSinTable[i]  = RoundFloatToInt((float)(scale*cos(2*M_PI*(float)i/numPoints)), numBits-1) << numBits;
    cosSinTable[i] |= RoundFloatToInt((float)(scale*sin(2*M_PI*(float)i/numPoints)), numBits-1) & mask;
  }

  return cosSinTable;
}


void DisposeFixedCosSinTable(PackedType* cosSinTable)
{
  if(cosSinTable)
    free(cosSinTable);
}

int GetFixedScalingBits(int N_FFT)
{
	int count = 0;
	while (N_FFT>2)
	{
		N_FFT = N_FFT>>1;
		count++;
	}
	return count;
}


int Radix4_CplxFixedFFT(const int* x, int* y,
                   const PackedType* cosSinTable,
                   int trigTableStrides, int N_FFT)
{
  int idx;
  int span;
  int twiddle;
  int strides;
	
  int x0_r, x0_i;
  int x1_r, x1_i;
  int x2_r, x2_i;
  int x3_r, x3_i;
  
  
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
    Radix4_Fixed_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);
    
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
      Radix4_Fixed_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);

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
        CplxMulConjFixed(x3_r, x3_i, x2_r, x2_i, cosSin3);
        CplxMulConjFixed(x2_r, x2_i, x1_r, x1_i, cosSin1);
        CplxMulConjFixed(x1_r, x1_i, x0_r, x0_i, cosSin2);
        
        x0_r = y[2*idx];
        x0_i = y[2*idx+1];
            
        /* W0 butterfly */
        Radix4_Fixed_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 1);
        
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
		Radix2_Fixed_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 1);
		
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
			
			CplxMulConjFixed(x1_r, x1_i, x0_r, x0_i, cosSin1);

			x0_r = y[2*idx];
			x0_i = y[2*idx+1];
						
			/* W0 radix 2 butterfly */
			Radix2_Fixed_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 0);
			
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
void CplxFixedToReal(int* x, const PackedType* cosSinTable, int N_FFT)
{
  int idx;
  
  int x0_r, x0_i;
  int x1_r, x1_i;
  int x2_r, x2_i;
  PackedType cosSin;
  
  /*
   X(0) = Re[X(0)] + Im[X(0)]
   X(N) = Re[X(0)] - Im[X(0)]
   Since X(0) and X(N) are real, we map X(N) to
   Im[X(N)]
   */
  x0_r = x[0];
  x0_i = x[1];
  
  x1_r = ScaleFxpRight(x0_r, 1) + ScaleFxpRight(x0_i, 1);
  x0_i = ScaleFxpRight(x0_r, 1) - ScaleFxpRight(x0_i, 1);
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
    
    CplxMulConjFixed(x2_r, x2_i, x1_r, x1_i, cosSin);
    
    x1_r = ScaleFxpRight(x0_r, 1) + x2_i;
    x1_i = ScaleFxpRight(x0_i, 1) - x2_r;
    x0_r = ScaleFxpRight(x0_r, 1) - x2_i;
    x0_i = -(ScaleFxpRight(x0_i, 1) + x2_r);
    
    x[2*idx]           = x1_r;
    x[2*idx+1]         = x1_i;
    x[2*(N_FFT-idx)]   = x0_r;
    x[2*(N_FFT-idx)+1] = x0_i;
    ++idx;
    
  } while(idx < N_FFT/2);
  
  /* Value at N_FFT/2: only change sign of imaginary part */
  x[2*idx]   = ScaleFxpRight(x[2*idx], 1);
  x[2*idx+1] = -ScaleFxpRight(x[2*idx+1], 1);

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
void ComplexFixedToRealInversion(const int* in, int* out, const PackedType* cosSinTable, int N_FFT)
{
	int idx, N_HALF;
	int x0_r, x0_i, x1_r, x1_i, x2_r, x2_i;
	int diag,onePlusSin,negOneMinSin,cosine;
	PackedType cosSin, cosVal, sinVal;

	N_HALF = N_FFT>>1;
	
	x0_r = in[0];
	x0_i = in[1];

	x1_r = ScaleFxpRight(x0_r, 0) + ScaleFxpRight(x0_i, 0);
	x0_i = ScaleFxpRight(x0_r, 0) - ScaleFxpRight(x0_i, 0);
	x0_r = x1_r;

	out[0] = ScaleFxpRight(x0_r, 0); // scale amount * 0.5
	out[1] = ScaleFxpRight(x0_i, 0);

	out[N_FFT] = in[N_FFT]<<1;
	out[N_FFT+1] = -(in[N_FFT+1]<<1);

	idx = 1;
  
  do
  {
	x0_r = in[2*idx];
	x0_i = in[2*idx+1];
	x1_r = in[2*(N_FFT-idx)];
	x1_i = in[2*(N_FFT-idx)+1];
	cosSin = cosSinTable[idx];
	// get everything into the most significant bits
	sinVal = (cosSin & 0x0000ffff)<<16;
	cosVal = (cosSin>>16)<<16;
	// diag = (1.0f - cosSin.y)*0.5f;
	diag = (0x40000000 - (sinVal>>1));
	// onePlusSin = (1.0f + cosSin.y)*0.5f;
	onePlusSin = (0x40000000 + (sinVal>>1));
	// negOneMinSin = (-1.0f - cosSin.y)*0.5f;
	negOneMinSin = (0xc0000000 - (sinVal>>1));
	// cosine
	cosine = cosVal>>1;

	// see float version below to see what the inverse matrix is
	// the up shift by two bits is to make up for the scaling in the
	// CplxToReal function.
	x2_r = FixedPointMul(diag,x0_r);
	x2_i = -FixedPointMul(cosine,x0_i);
	x2_r = FixedPointMulAndAdd(onePlusSin, x1_r, x2_r);
	x2_i = FixedPointNegMulAndAdd(cosine,x1_i,x2_i);
	out[2*idx] = (x2_r + x2_i)<<2;
	x2_r = FixedPointMul(onePlusSin,x0_r);
	x2_r = FixedPointMulAndAdd(diag,x1_r,x2_r);
	out[2*(N_FFT-idx)] = (x2_r - x2_i)<<2;
	x2_r = FixedPointMul(cosine,x0_r);
	x2_i = FixedPointMul(diag,x0_i);
	x2_r = FixedPointNegMulAndAdd(cosine,x1_r,x2_r);
	x2_i = FixedPointMulAndAdd(negOneMinSin,x1_i,x2_i);
	out[2*idx+1] = (x2_r + x2_i)<<2;
	x2_i = FixedPointMul(negOneMinSin,x0_i);
	x2_i = FixedPointMulAndAdd(diag,x1_i,x2_i);
	out[2*(N_FFT-idx)+1] = (x2_r + x2_i)<<2;

    ++idx;
  } while(idx < N_HALF);

}

void InvertFixedOutput(int* x, int N_FFT)
{
	int i;
	float scale = 1.0f / N_FFT;
	
	/* the first pair are merely scaled */
	x[0] = (int)((float)x[0] * scale);
	x[1] = (int)((float)x[1] * scale);
	
	/* the rest are scaled and swapped about the middle */
	for (i = 1; i <= N_FFT / 2; i++)
	{
		int a_index = 2 * i;
		int b_index = 2 * (N_FFT - i);
		
		int a_real = x[a_index];
		int a_imag = x[a_index + 1];
		int b_real = x[b_index];
		int b_imag = x[b_index + 1];
		
		x[a_index]     = (int)(b_real * scale);
		x[a_index + 1] = (int)(b_imag * scale);
		x[b_index]     = (int)(a_real * scale);
		x[b_index + 1] = (int)(a_imag * scale);
	}
}

