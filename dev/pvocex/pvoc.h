/* Functions in PVOC.C */

void    warpse(float*,float*,int,double);
void    usage(void);
void    hamming(float*,int,int);
float   *float_array(int);
void    malerr(char*,int);
void    kaiser_(int*,float*,int*,int*,float*);


/* Functions in MXFFT.C */

void fft_(float *, float *,int,int,int,int);
void fftmx(float *,float *,int,int,int,int,int,int *,float *,float *,float *,float *,int *,int[]);
void reals_(float *,float *,int,int);

#ifndef min
#define min(a,b)    (((a)<(b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a)>(b)) ? (a) : (b))
#endif
