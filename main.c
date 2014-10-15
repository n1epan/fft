#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <io.h>
#define snprintf _snprintf
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#endif

#include "fft.h"

#define MACRO_TO_STR(s) MACRO_TO_STR2(s)
#define MACRO_TO_STR2(s) #s

typedef int (*TestFunction)(void);

/* tests */
static int test_all(void);
static int eric_test(void);
static int inverse_test(void);

static struct
{
	const char* arg;
	const char* name;
	TestFunction test;
}
g_test_table[] =
{
	{ "-test_all", "Run all tests", test_all },
	{ "-eric", "Eric Allamanche test", eric_test },
	{ "-inverse", "Inverse FFT test", inverse_test }
};

static int g_test_count = sizeof(g_test_table) / sizeof(g_test_table[0]);
static int g_on_console = 0;

static void print_usage();
static void run_arg_test(const char* arg);
static void run_test(int test_index);

int main(int argc, char** argv)
{
	int i;
	
	/* initialize */
	g_on_console = isatty(fileno(stdout));
	
	printf(
		"FFT Library Test program (%d bit, using %s data type)\n",
		(int)(8 * sizeof(void*)),
		MACRO_TO_STR(Type));
	
	/* handle usage display */
	if (argc < 2)
	{
		print_usage();
		
		return 0;
	}
	
	/* execute arguments */
	for (i = 1; i < argc; i++)
	{
		run_arg_test(argv[i]);
	}
	
	return 0;
}

static void print_usage()
{
	int i;
	
	printf("usage: fft_test <test_name>\n");
	
	printf("Tests:\n");
	
	for (i = 0; i < g_test_count; i++)
	{
		printf(
			"  %-14s  %s\n",
			g_test_table[i].arg,
			g_test_table[i].name);
	}
}

static void run_arg_test(const char* arg)
{
	int i;
	
	for (i = 0; i < g_test_count; i++)
	{
		if (strcmp(arg, g_test_table[i].arg) == 0)
		{
			run_test(i);
		}
	}
}

static void run_test(int test_index)
{
	int result;
	
	printf("%s . . .\n", g_test_table[test_index].name);
	
	result = g_test_table[test_index].test();
	
	printf(
		"%s: %s\n",
		g_test_table[test_index].name,
		result ? "Pass" : "FAIL");
}

static int test_all(void)
{
	int i;
	
	for (i = 0; i < g_test_count; i++)
	{
		if (g_test_table[i].test != test_all)
		{
			run_test(i);
		}
	}
	
	return 1;
}

static int eric_test(void)
{
	const int M = 15;
	const int N_FFT = 1<<(M);

	double maxErr = 0.0;
	float scale;
	float sign;
	double K_r = 0.0;
	double K_i = 0.0;
	double e1, e2;

	Type* x = (Type*)malloc(sizeof(Type)*2*N_FFT);
	Type* y = (Type*)malloc(sizeof(Type)*2*N_FFT);

	scale = 1.0; //pow(2, -31);

	printf("FFT size: %d\n", 2*N_FFT);

#if 0
	{
		int i;

		HCplxFFT hFFT;
		
		CreateCplxFFT(&hFFT, N_FFT);

		srand(0);
		for(i = 0; i < N_FFT; ++i)
		{
			//x[2*i] = 0; //Round(cos(2*M_PI*1*i/N_FFT), 31);
			//x[2*i+1] = 0;
			x[2*i] = (rand()<<1)>>1;
			x[2*i+1] = (rand()<<1)>>1;
		}
		//x[2*1] = (1<<31)-1; //Round(1, 31);

		ComputeCplxFFT(hFFT, x, y, Forward);

		{
			printf("\nReference FFT:\n");
			int k;
			int n;
			double maxErr = 0.0;
			for(k = 0; k < N_FFT; ++k)
			{
				double K_r = 0.0;
				double K_i = 0.0;
				double e1, e2;
				for(n = 0; n < N_FFT; ++n)
				{
					K_r += scale*x[2*n]*cos(2*M_PI*n*k/N_FFT)+scale*x[2*n+1]*sin(2*M_PI*n*k/N_FFT);
					K_i += scale*x[2*n+1]*cos(2*M_PI*n*k/N_FFT)-scale*x[2*n]*sin(2*M_PI*n*k/N_FFT);
				}
				K_r /= N_FFT;
				K_i /= N_FFT;
				e1 = K_r*K_r + K_i*K_i;
				e2 = (scale*y[2*k])*(scale*y[2*k]) + (scale*y[2*k+1])*(scale*y[2*k+1]);
				//printf("x_fix: %5.4f\tx_ref: %5.4f\ty_fix: %5.4f\ty_ref: %5.4f\t (err: %5.4g)\n", scale*y[2*k], K_r, scale*y[2*k+1], K_i, fabs(e1-e2));
				if(fabs((e1-e2)/e1) > maxErr)
					maxErr = fabs((e1-e2)/e1);

				//printf("i: %2d\ty_r: %f\ty_i: %f\n", k, fabs(K_r-y[2*k]), fabs(K_i-y[2*k+1]));
				//printf("i: %2d\ty_r: %f\ty_i: %f\n", k, K_r, K_i);
			}
			printf("Max relative error: %6.5g\n", maxErr);
		}  

		DisposeCplxFFT(hFFT);
	}
#endif

	{
		int i,n,k;

		/* Real FFT */
		HCplxFFT hFFT;

		printf("Creating signal...\n");

		//srand(0);
		for(i = 0; i < 2*N_FFT; ++i)
		{
			x[i] = (Type)(rand()<<1)/pow(2,31);
			//x[i] = 0X7fffffffL + (rand() & 0x01L);
		}
		//x[0] = 1.0f; //(1<<30);

		printf("Computing real-valued FFT...\n");
		CreateRealFFT(&hFFT, 2*N_FFT);

		ComputeRealFFT(hFFT, x, y, Forward);

		DisposeRealFFT(hFFT);

		//for(i = 0; i < N_FFT; ++i)
		//printf("x: %7.4g\ty: %7.4g\n", scale*y[2*i], scale*y[2*i+1]);

		printf("Comparing with reference FT (not fast)...\n");

		/* Reference FFT */
		sign = 1;
		K_r = K_i = 0.0f;
		for(n = 0; n < 2*N_FFT; ++n)
		{
			K_r += scale*x[n];
			K_i += sign*scale*x[n];
			sign = -sign;
		}
		//K_r /= 2*N_FFT;
		//K_i /= 2*N_FFT;
		//printf("x: %7.4g\ty: %7.4g\n", K_r, K_i);
		e1 = K_r*K_r + K_i*K_i;
		e2 = (scale*y[0])*(scale*y[0]) + (scale*y[1])*(scale*y[1]);

		if(fabs((e1-e2)) > maxErr)
			maxErr = fabs((e1-e2));

		for(k = 1; k < N_FFT; ++k)
		{
			if (g_on_console)
			{
				printf("\r%.1f%%", 100.0 * k / N_FFT);
				fflush(stdout);
			}
			
			K_r = K_i = 0.0f;
			for(n = 0; n < 2*N_FFT; ++n)
			{
				K_r += scale*x[n]*cos(2*M_PI*n*k/(2*N_FFT));
				K_i += -scale*x[n]*sin(2*M_PI*n*k/(2*N_FFT));
			}
			//K_r /= 2*N_FFT;
			//K_i /= 2*N_FFT;
			//printf("x: %7.4g\ty: %7.4g\n", K_r, K_i);
			e1 = K_r*K_r + K_i*K_i;
			e2 = (scale*y[2*k])*(scale*y[2*k]) + (scale*y[2*k+1])*(scale*y[2*k+1]);

			if(fabs((e1-e2)) > maxErr)
				maxErr = fabs((e1-e2));
		}
		
		if (g_on_console)
		{
			printf("\n");
		}

		printf("Max relative error: %8.7e\n", maxErr);
	}

	free(y);  
	free(x);

	return 1;
}

static int inverse_test(void)
{
	const int size = 16;
	Type x[2 * size] = { 0 };
	Type y[2 * size] = { 0 };
	Type z[2 * size] = { 0 };
	HCplxFFT h = 0;
	int i;
	double sum = 0;
	double rms = 0;
	
	/* create a test signal */
	for (i = 0; i < size; i++)
	{
		x[2 * i]     = (Type)(rand()<<1)/pow(2,31);
		x[2 * i + 1] = (Type)(rand()<<1)/pow(2,31);
		//x[2 * i]     = (Type)sin(2 * M_PI * i / size);
		//x[2 * i + 1] = (Type)cos(2 * M_PI * i / size);
	}
	
	/* create the FFT object and perform forward and inverse FFTs */
	CreateCplxFFT(&h, size);

	ComputeCplxFFT(h, x, y, Forward);
	
	ComputeCplxFFT(h, y, z, Inverse);

	/* display the results */
	for (i = 0; i < size; i++)
	{
		Type diff_real = x[2 * i] - z[2 * i];
		Type diff_imag = x[2 * i + 1] - z[2 * i + 1];
		
		sum += diff_real * diff_real + diff_imag * diff_imag;
		
		printf(
			"x[%2d] = (%6.3f, %6.3f);  y[%2d] = (%6.3f, %6.3f);  z[%2d] = (%6.3f, %6.3f)\n",
			i,
			(double)x[2 * i],
			(double)x[2 * i + 1],
			i,
			(double)y[2 * i],
			(double)y[2 * i + 1],
			i,
			(double)z[2 * i],
			(double)z[2 * i + 1]);
	}
	
	rms = sqrt(sum / (size * 2));
	
	printf("RMS difference: %G\n", rms);
	
	/* clean up */
	DisposeCplxFFT(h);
	
	return (rms < 1e-3);
}

