/*
 * Copyright (C) 2005 Gracenote, Inc. All Rights Reserved.
 *     
 */


#include <stdlib.h>
#include <math.h>
#include "fft_impl.h"
#include "functions.h"
#ifdef FXP_32_BIT
#include "fft_bfly_fix.h"
#else
#include "fft_bfly_float.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

FILE* g_out = NULL;


PackedType* CreateCosSinTable(int numPoints, int numBits)
{
#ifdef FXP_32_BIT
  const int mask = (1<<(numBits))-1;
  float scale = (float)pow(2.0, numBits-1);
#endif
  int i;
  
  int size = 3*(numPoints/4 - 1) + 1;
  
  PackedType* cosSinTable = (PackedType*)malloc(sizeof(PackedType)*size);
  
  if(NULL == cosSinTable)
    return NULL;
  
  for(i = 0; i < size; ++i)
  {
#ifdef FXP_32_BIT
	//HeeYoung 5/12/06: casted to float to remove warning
    cosSinTable[i]  = RoundFloatToInt((float)(scale*cos(2*M_PI*(float)i/numPoints)), numBits-1) << numBits;
    cosSinTable[i] |= RoundFloatToInt((float)(scale*sin(2*M_PI*(float)i/numPoints)), numBits-1) & mask;
#else
    cosSinTable[i].x = cos(2*M_PI*(float)i/numPoints);
    cosSinTable[i].y = sin(2*M_PI*(float)i/numPoints);
#endif
  }

  return cosSinTable;
}


void DisposeCosSinTable(PackedType* cosSinTable)
{
  if(cosSinTable)
    free(cosSinTable);
}

int GetScalingBits(int N_FFT)
{
	int count = 0;
	while (N_FFT>2)
	{
		N_FFT = N_FFT>>1;
		count++;
	}
	return count;
}


int Radix4_CplxFFT(const Type* x, Type* y,
                   const PackedType* cosSinTable,
                   int trigTableStrides, int N_FFT)
{
  int idx;
  int span;
  int twiddle;
  int strides;
	
  Type x0_r, x0_i;
  Type x1_r, x1_i;
  Type x2_r, x2_i;
  Type x3_r, x3_i;
  
  
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
    Radix4_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);
    
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
      Radix4_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 2);

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
        CplxMulConj(x3_r, x3_i, x2_r, x2_i, cosSin3);
        CplxMulConj(x2_r, x2_i, x1_r, x1_i, cosSin1);
        CplxMulConj(x1_r, x1_i, x0_r, x0_i, cosSin2);
        
        x0_r = y[2*idx];
        x0_i = y[2*idx+1];
            
        /* W0 butterfly */
        Radix4_Bfly_W0(x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i, 1);
        
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
		Radix2_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 1);
		
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
			
			CplxMulConj(x1_r, x1_i, x0_r, x0_i, cosSin1);

			x0_r = y[2*idx];
			x0_i = y[2*idx+1];
						
			/* W0 radix 2 butterfly */
			Radix2_Bfly_W0(x0_r, x0_i, x1_r, x1_i, 0);
			
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
void CplxToReal(Type* x, const PackedType* cosSinTable, int N_FFT)
{
  int idx;
  
  Type x0_r, x0_i;
  Type x1_r, x1_i;
  Type x2_r, x2_i;
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
    
    CplxMulConj(x2_r, x2_i, x1_r, x1_i, cosSin);
    
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
void ComplexToRealInversion(const Type* in, Type* out, const PackedType* cosSinTable, int N_FFT)
{
	int idx, N_HALF;
	Type x0_r, x0_i, x1_r, x1_i, x2_r, x2_i;
	Type diag,onePlusSin,negOneMinSin,cosine;
	PackedType cosSin, cosVal, sinVal;
	Type scale;

	N_HALF = N_FFT>>1;
	
#ifdef FXP_32_BIT
	x0_r = in[0];
	x0_i = in[1];

	x1_r = ScaleFxpRight(x0_r, 0) + ScaleFxpRight(x0_i, 0);
	x0_i = ScaleFxpRight(x0_r, 0) - ScaleFxpRight(x0_i, 0);
	x0_r = x1_r;

	out[0] = ScaleFxpRight(x0_r, 0); // scale amount * 0.5
	out[1] = ScaleFxpRight(x0_i, 0);

	out[N_FFT] = in[N_FFT]<<1;
	out[N_FFT+1] = -(in[N_FFT+1]<<1);
#else
	scale = 1.0f/(float)2.0f*N_FFT;
	out[0] = (in[0] + in[1])*scale*0.5f;
	out[1] = (in[0] - in[1])*scale*0.5f;
	out[N_FFT] = in[N_FFT]*scale;
	out[N_FFT+1] = -in[N_FFT+1]*scale;
#endif

	idx = 1;
  
  do
  {
	x0_r = in[2*idx];
	x0_i = in[2*idx+1];
	x1_r = in[2*(N_FFT-idx)];
	x1_i = in[2*(N_FFT-idx)+1];
#ifdef FXP_32_BIT
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
#else
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
#endif
    ++idx;
  } while(idx < N_HALF);
}

/*
// Split radix Inverse FFT - Sorenson, et al - not tested, but possible faster implementation - BobC
void RealInverseFFT(const Type* input, Type* out, const PackedType* cosSinTable, int N_FFT)
{
	int i,j,k,i5,i6,i7,i8,i0,id,i1,i2,i3,i4,N_DOUBLE,N_4TH,N_8TH,N_MINUS1,N_HALF;
	float tmp1,tmp2,tmp3,tmp4,tmp5,a3,sin1,sin3,cos1,cos3,a,e,root2,scale;

	root2 = sqrtf(2.0f);
	scale = 0.0f;//1.0f/N_FFT;

	N_MINUS1 = N_FFT-1;
	N_DOUBLE = N_FFT<<1; // double FFT size
	N_HALF = N_FFT>>1; // half FFT size
	

	// Step 1/2 Get the input from real(1),imaj(1),real(2),imaj(2),...
	// to real(1), real(2),...real(N/2-1),imaj(N/2-1), ... imaj(2), imaj(1)
	// move the evens
	for (i=0; i<N_HALF; i=i++)
	{
		out[i] = input[i*2]*scale;
		out[i+N_HALF] = input[N_MINUS1-i*2]*scale;
	}
	
	// Step 1 Do the inverse FFT
	for (k=N_FFT; k>2; k>>=1)
	{
		id=N_DOUBLE;
		N_DOUBLE>>=1; // divide by 2
		N_4TH=N_DOUBLE>>2; //divide by 4
		N_8TH=N_DOUBLE>>3; //divide by 8
		e = 2.0f*(float)M_PI/(float)N_DOUBLE;
		i1 = 0;

		do 
		{
			// i1 takes the value set at the bottom of the do statement
			for (; i1<N_FFT; i1+=id)
			{
				i2 = i1 + N_4TH;
				i3 = i2 + N_4TH;
				i4 = i3 + N_4TH;
				tmp1 = out[i1] - out[i3];
				out[i1] += out[i3];
				out[i2] *= 2.0f;
				out[i3] = tmp1 - 2.0f*out[i4];
				out[i4] = tmp1 + 2.0f*out[i4];
				if (N_4TH != 1)
				{
					i0 = i1 + N_8TH;
					i2 += N_8TH;
					i3 += N_8TH;
					i4 += N_8TH;
					tmp1 = (out[i2] - out[i0]) / root2;
					tmp2 = (out[i3] + out[i4]) / root2;
					out[i0] += out[i2];
					out[i2] = out[i4] - out[i3];
					out[i3] = 2.0f * (-tmp1 - tmp2);
					out[i4] = 2.0f * (tmp1 - tmp2);
				}
			}
			id <<= 1; // mult by 2
			i1 = id - N_DOUBLE;
			id <<= 1;
		} while (i1 < N_MINUS1);
		a = e;
		for (j=2; j <= N_8TH; j++)
		{
			a3 = 3 * a;
			cos1 = cosf(a);
			sin1 = sinf(a);
			cos3 = cosf(a3);
			sin3 = sinf(a3);
			a = j*e;
			i = 0;
			id = N_DOUBLE<<1;
			do
			{
				for (; i<N_FFT; i += id)
				{
					i1 = i + j -1;
					i2 = i1 + N_4TH;
					i3 = i2 + N_4TH;
					i4 = i3 + N_4TH;
					i5 = i4 + N_4TH;
					i6 = i5 + N_4TH;
					i7 = i6 + N_4TH;
					i8 = i7 + N_4TH;
					tmp1 = (out[i1] - out[i6]);
					out[i1] += out[i6];
					tmp2 = (out[i5] - out[i2]);
					out[i5] += out[i2];
					tmp3 = out[i8] + out[i3];
					out[i6] = out[i8] - out[i3];
					tmp4 = out[i4] + out[i7];
					out[i2] = out[i4] - out[i7];
					tmp5 = tmp1 - tmp4;
					tmp1 += tmp4;
					tmp4 = tmp2 - tmp3;
					tmp2 += tmp3;
					out[i3] = tmp5 * cos1 + tmp4 * sin1;
					out[i7] = tmp4 * cos1 + tmp5 * sin1;
					out[i4] = tmp1 * cos3 - tmp2 * sin3;
					out[i8] = tmp2 * cos3 + tmp1 * sin3;
				} // for i
				id <<= 1; // mult by 2
				i = id - N_DOUBLE;
				id <<= 1;
			} while(i < N_MINUS1);
		} // for j = 2
	} // for k=N_FFT;

	// Step 2 - Radix 2 butterfly
	i0 = 0;
	id = 4;
	do 
	{
		for (; i0<N_MINUS1; i0+=id)
		{
			i1 = i0 + 1;
			tmp1 = out[i0];
			out[i0] = tmp1 + out[i1];
			out[i1] = tmp1 - out[i1];
		}
		id <<= 1; // mult by 2
		i0 = id - 2;
		id <<= 1;
	} while (i0 < N_MINUS1);

	// Step 3: Get the output back to normal order
	for (i=0, j=0, N_DOUBLE=N_FFT/2; i<N_MINUS1; i++)
	{
		if (i < j)
		{
			tmp1 = out[j];
			out[j] = out[i];
			out[i] = tmp1;
		}
		k = N_DOUBLE;
		while (k <= j)
		{
			j -= k;
			k >>= 1; // divide by 2
		}
		j += k;
	} // for i=0; etc.
}

void RotateSpectrum (Type* x, int N_FFT)
{
	int i,N_HALF,aIdx,bIdx;
	Type aReal,aImaj;
	N_HALF = N_FFT>>1;

	for (i=1; i<N_HALF; i++)
	{
		aIdx = 2*i;
		bIdx = 2*(N_FFT-i);
		aReal = x[aIdx];
		aImaj = x[aIdx+1];
		
		x[aIdx]		= x[bIdx];
		x[aIdx+1]	= x[bIdx+1];
		x[bIdx]		= aReal;
		x[bIdx+1]	= aImaj;
	}
}
*/

