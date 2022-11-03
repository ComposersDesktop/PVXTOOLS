/* phase vocoder */
/* Adapted by Richard Dobson to use and test the new pvocex file format.
 */

/* -----------------------------------------------------------------------------
    HISTORY.
    This code is based on pvoc.c from the CARL distribution, as adapted for use in the
    CDP system. It supports only mono files.
    There is one inherited problem: 
    because the analysis operation puts a rounded number of analysis
    samples into the analysis file, the EXACT original file length
    cannot be restored in the newly synthesised file. While this
    may be annoying, and could no doubt be corrected, it is not
    likely to be of critical importance, since the difference is
    less than one analysis window in samples.

    The fft was 'translated' from fortran to C by Trevor Wishart and Andrew Bentley.
    The fft was made more efficient by Keith Henderson.
    
    the code also incorporates CDP contributions from Martin Atkins and Nick Laviers
    

------------------------------------------------------------------------------*/

/* PVOCEX version history
 * May 2000: v0.1 initial version.
 * June 18 2000 incorporated M (analWinlen) correctly.
 * v0.12:
 * added correction to sinc function code from Dan Timis 
 * July 15 2000  added fftw support for huge speed increase, at cost of larger binary
 * July 20 2000  added -m flag for all those bad unix audio apps,
 *               added -v flag - progress messages do slow things down a little.
 * v0.13
 * march 2001    added -B and -x flags (mainly for use with impulse responses)
 *               removed limit on fft sizes (for the same purpose)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef unix
#include <sys/timeb.h>
#include <unistd.h>
#else
#include <sys/timeb.h>
#endif
#include "pvoc.h"   
#include "pvfileio.h"
#include "cmdutils.h"
#include <portsf.h>

/*RWD Oct12 */
#include <assert.h>

# ifdef USE_FFTW
# include <rfftw.h>
# endif

double
timer();
void            
stopwatch(int flag);
#define PROG    "PVOCEX"

#define LOSR    (44100)     
#undef TWOPI            /* Goodness knows why ??? */

#ifndef PI
#define PI  (3.14159265358979323846)
#endif

#define TWOPI   (2*PI)
/*RWD 3:2000 this rapid report slows things down a lot!: was 0.01 */
#define TIME_INTERVAL   0.5 /*Time between progress report messages*/
#define MONO    (1)     /*In absence of defn. in sfsys.h*/

#ifndef  WAVE_FORMAT_IEEE_FLOAT
#define WAVE_FORMAT_IEEE_FLOAT (0x0003)
#endif



#ifndef NOOVERCHK
int num_overflows = 0;
double maxsample = 0.0;
double minsample = 0.0;
#endif

void help(void);
/* use for profiling only! */
static float myhypot(float x,float y)
{
    return (float) sqrt((double)(x*x  + y*y));
}

static float myatan2(double i, double r)
{
    return (float) atan2(i,r);
}


void
warpse(float *anal,float *env,int N2,double warp)
{/* spectral envelope detection: this is a very crude peak picking algorithm
    which is used to detect and pre-warp the spectral envelope so that
    pitch transposition can be performed without altering timbre.
    The basic idea is to disallow large negative slopes between
    successive values of magnitude vs. frequency. */

    float   lastmag, mag, nextmag, slope, pkOld,
        eps = (float)(-32.0/N2);
    int pkcnt, i, j;

    lastmag = *anal;
    mag = *(anal + 2);
    pkOld = lastmag;
    *env = pkOld;
    pkcnt = 1;

    for (i = 1; i <= N2; i++){  /* step thru spectrum*/

        nextmag = (float)( i<N2 ? anal[2*i+2] : 0.0 );

        slope = (float)(pkOld != 0.0 ? (mag-pkOld)/(pkOld*pkcnt) : -10.0 );

                        /* look for peaks */
        if ((mag>=lastmag)&&(mag>nextmag)&&(slope>eps)){
            *(env + i) = mag;
            pkcnt--;
            for (j = 1; j <= pkcnt; j++){
                *(env + i - pkcnt + j - 1)
                = (float)(pkOld * (1. + slope * j));
            }
            pkOld = mag;
            pkcnt = 1;
        }
        else pkcnt++;       /* not a peak */

        lastmag = mag;
        mag = nextmag;
    } /* end of scan through spectrum */

    if (pkcnt > 1) {        /* get final peak */
        mag = *(anal + N2*2);
        slope = ((float) (mag - pkOld) / pkcnt);
        *(env + N2) = mag;
        pkcnt--;
        for (j = 1; j <= pkcnt; j++){
            *(env + N2 - pkcnt + j - 1) = pkOld + slope * j;
        }
    }

    for (i = 0; i <= N2; i++){  /* warp spectral env.*/
        j = (int)((float) i * warp);
        if ((j <= N2) && (*(env + i) != 0.))
        *(anal + 2*i) *= *(env + j) / *(env + i);
        else *(anal + 2*i) = 0.0f;
    }
}

void
usage(void)
{
        fprintf(stderr,"USAGE: pvocex [flags] infile outfile\n\n");
        fprintf(stderr,"\ttype pvocex -h for more information on flags.\n");
        exit(1);
}

void help(void)
{
        fprintf(stderr,"R = input sample rate (automatically read from stdin)\n");
        fprintf(stderr,"F = fundamental frequency (R/256) DONT'T USE -F AND -N\n");
        fprintf(stderr,"N = sample window length for fft (256 unless -F is specified)\n");
        fprintf(stderr,"W = filter overlap factor: {0,(1),2,3} DON'T USE -W AND -M\n");
        fprintf(stderr,"M = analysis window length (N*2 unless -W is specified)\n");
        fprintf(stderr,"L = synthesis window length (M) \n");
        fprintf(stderr,"D = decimation factor (min((M/(8*T)),(M/8))\n");
        fprintf(stderr,"I = interpolation factor (=T*D) \n");
        fprintf(stderr,"T = time-scale factor (1.)\n");
        fprintf(stderr,"P = pitch-scale factor (1.) DON'T USE -T AND -P\n");
        fprintf(stderr,"C = resynthesize odd (1) or even (2) channels only\n");
        fprintf(stderr,"i = resynthesize bandpass filters i thru j only\n");
        fprintf(stderr,"j = resynthesize bandpass filters i thru j only\n");
        fprintf(stderr,"b = starting sample (0)\n");
        fprintf(stderr,"e = final sample (end of input)\n");
        fprintf(stderr,"w = warp factor for spectral envelope (1.)\n");
        fprintf(stderr,"A:  analysis only: output analysis data (.pvx format)\n");
        fprintf(stderr,"x:  write data as PVOC_COMPLEX (default: AMP_FREQ format)\n");
#ifdef NOTDEF
        /*RWD these options not supported yet */
        fprintf(stderr,"E:  analysis only: output spectral envelope\n");
        fprintf(stderr,"X:  analysis only: output magnitude values\n");
#endif
        fprintf(stderr,"S:  synthesis only: input analysis data (.pvx format)\n");
        fprintf(stderr,"K:  use Kaiser Window \n");
        fprintf(stderr,"B:  use Box (rectangular) Window\n");
        fprintf(stderr,"    (default Window: Hamming)\n");
        fprintf(stderr,"m:  write minimal file header (no PEAK data)\n");
        fprintf(stderr,"v:  suppress progress messages\n");
#ifdef NOTDEF
        /*RWD hidden for now */
        fprintf(stderr,"V [filename]:  verbose (summarize on pvoc.s or file)\n");
        fprintf(stderr,"if filename is specified, it must follow all flags\n");
#endif

#ifdef _DEBUG
        fprintf(stderr,"DEBUG only: -pN:  generate text dumps of key arrays at sample N\n");
        fprintf(stderr,"            -H[-]N: dump  phase of bin N to binNphase.dat\n");
        fprintf(stderr,"            -Q[-]N: dump  freq of bin N to binNfreq .dat\n");
        fprintf(stderr,"                 (negative N indicates output, otherwise input.)\n");
#endif


        exit(1);
}



FUNCTION *gettsf(char *filename)
{   /* Gets pairs of values from the named file and sets up
       pointer to function header */

    int i=0,len,j,k,l,toggle = 1;
    double sample, *f, *x, *y;
    FUNCTION *p;
    FILE *fp;
    
    /* Open named file for time stretching data input */
    if( (fp = fopen(filename,"r")) == NULL){
        fprintf(stderr,"pvocex: Cannot open data file\n");
        return (NULL);
    }
    
    /* Read stream of values into buffer */
    while(fscanf(fp,"%lf",&sample) > 0 ){
        if(i == 0){
            f = (double *) malloc(sizeof(double)*BUFSIZ);
            if(f == NULL){
                fprintf(stderr,
                    "pvocex: cannot get memory for buffer\n");
                return (NULL);
            }
            len = BUFSIZ;
        }
        if(i < len-1)
            f[i++]= sample;
        else{
            len += BUFSIZ;
            f = (double *) realloc(f,len*sizeof(double));
            if(f == NULL){
                fprintf(stderr,
                    "pvocex: cannot get memory for buffer\n");
                return (NULL);
            }
            f[i++] = sample;
        }
    }

    /* Rearrange values into arrays */
    f = (double *) realloc(f, i*sizeof(double));
    if(f == NULL){
        fprintf(stderr,
            "pvocex: cannot get memory for buffer\n");
        return (NULL);
    }
    p = (FUNCTION *) calloc(sizeof(FUNCTION),1);
    if(p == NULL){
        fprintf(stderr,
            "pvocex: cannot get memory for buffer\n");
        return (NULL);
    }
    if((i%2) != 0){
        fprintf(stderr,"pvocex: time stretching data not in pairs\n");
        return(NULL);
    }
    x = (double *) malloc(sizeof(double)*i/2);
    y = (double *) malloc(sizeof(double)*i/2);
    for(j=k=l=0; j<i; j++,toggle = 1-toggle){
        if(toggle)
            x[k++] = f[j];
        else
            y[l++] = f[j];
    }
    p->fxval = x;
    p->fyval = y;
    p->flen = i/2;
    p->ftype = "XY pairs";
    p->fname = filename;

    free(f);
    return (p); 
}

void recwin(float *win,int winLen,int even)
{
    int i;
    if(even){
        for(i=0;i < winLen;i++)
            win[i] = 1.0f;
    }
    else{
        for(i=0;i <= winLen;i++)
            win[i] = 1.0f;
    }
}


void
hamming(float *win,int winLen,int even)
{
    float Pi,ftmp;
    int i;

/***********************************************************
                    Pi = (float)((double)4.*atan((double)1.));
***********************************************************/
    Pi = (float)PI;
    ftmp = Pi/winLen;

    if (even) {
        for (i=0; i<winLen; i++)
        *(win+i) = (float)((double).54 + (double).46*cos((double)(ftmp*((float)i+.5))));
//      *(win+winLen) = 0.0f;
    }
    else{   *(win) = 1.0f;
        for (i=1; i<=winLen; i++)
        *(win+i) =(float)((double).54 + (double).46*cos((double)(ftmp*(float)i)));
    }
    return;
}

void
vonhann(float *win,int winLen,int even)
{
    float Pi,ftmp;
    int i;

/***********************************************************
                    Pi = (4.*atan(1.));
***********************************************************/
    Pi = (float)PI;
    ftmp = Pi/winLen;

    if (even) {
        for (i=0; i<winLen; i++)
        *(win+i) = (float)(.5 + .5 *cos(ftmp*((double)i+.5)));
//      *(win+winLen) = 0.0f;
    }
    else{   *(win) = 1.0f;
        for (i=1; i<winLen; i++)  /*RWD Nov29 06 was <= */
        *(win+i) =(float)(.5 + .5 *cos(ftmp*(double)i));
    }
    return;
}

/******************** KAISER WINDOW ***************/
double besseli( double x)
{
    double ax, ans;
    double y;

    if (( ax = fabs( x)) < 3.75)     {
    y = x / 3.75;
    y *= y;
    ans = (1.0 + y * ( 3.5156229 +
              y * ( 3.0899424 +
                y * ( 1.2067492 +
                      y * ( 0.2659732 +
                        y * ( 0.360768e-1 +
                          y * 0.45813e-2))))));
    }
    else {
    y = 3.75 / ax;
    ans = ((exp ( ax) / sqrt(ax))
        * (0.39894228 +
           y * (0.1328592e-1 +
            y * (0.225319e-2 +
             y * (-0.157565e-2 +
                  y * (0.916281e-2 +
                   y * (-0.2057706e-1 +
                    y * (0.2635537e-1 +
                         y * (-0.1647633e-1 +
                          y * 0.392377e-2)))))))));
    }
    return ans;
}

//courtesy of Csound....

void kaiser(float *win,int len,double beta)
{
    float *ft = win;
    double i,xarg = 1.0;    //xarg = amp scalefactor
    for (i = -len/2.0 + .1 ; i < len/2.0 ; i++)
        *ft++ = (float) (xarg *
              besseli((beta * sqrt(1.0-pow((2.0*i/(len - 1)), 2.0)))) /
              besseli( beta));
   // assymetrical hack: sort out first value!
   win[0] = win[len-1];
}

float
*float_array(int nnn)
{   /* set up a floating point array length nnn. */
    float *ptr;
    ptr = (float *) calloc(nnn,sizeof(float));
    if(ptr==NULL){
        fprintf(stderr,"pvocex: insufficient memory\n");
        exit(1);
    }
    return (ptr);
}

void
malerr(char *str, int ex)
{
    fprintf(stderr, "%s\n", str);
    exit(1);
}




int
main(int argc, char *argv[])
{
    double  rratio;
    float   *input,     /* pointer to start of input buffer */
        *output,    /* pointer to start of output buffer */
        *anal,      /* pointer to start of analysis buffer */
        *syn,       /* pointer to start of synthesis buffer */
        *banal,     /* pointer to anal[1] (for FFT calls) */
        *bsyn,      /* pointer to syn[1]  (for FFT calls) */
        *nextIn,    /* pointer to next empty word in input */
        *nextOut,   /* pointer to next empty word in output */
        *analWindow,    /* pointer to center of analysis window */
        *synWindow, /* pointer to center of synthesis window */
        *maxAmp,    /* pointer to start of max amp buffer */
        *avgAmp,    /* pointer to start of avg amp buffer */
        *avgFrq,    /* pointer to start of avg frq buffer */
        *env,       /* pointer to start of spectral envelope */
        *i0,        /* pointer to amplitude channels */
        *i1,        /* pointer to frequency channels */
        *oi,        /* pointer to old phase channels */
        *oldInPhase,    /* pointer to start of input phase buffer */
        *oldOutPhase;   /* pointer to start of output phase buffer */

    int m, n;

    int N = 0,      /* number of phase vocoder channels (bands) */
        M = 0,      /* length of analWindow impulse response */
        L = 0,      /* length of synWindow impulse response */
        D = 0,      /* decimation factor (default will be M/8) */
        I = 0,      /* interpolation factor (default will be I=D)*/
        W = -1,     /* filter overlap factor (determines M, L) */
    /*  F = 0,  */  /* fundamental frequency (determines N) */
    /*  F2 = 0, */  /* F/2 */
        analWinLen, /* half-length of analysis window */
        synWinLen,  /* half-length of synthesis window */       
        outCount,   /* number of samples written to output */
        ibuflen,    /* length of input buffer */
        obuflen,    /* length of output buffer */
        nI = 0,     /* current input (analysis) sample */
        nO,     /* current output (synthesis) sample */
        nMaxOut,    /* last output (synthesis) sample */
        nMin = 0,   /* first input (analysis) sample */
        nMax = 100000000;   /* last input sample (unless EOF) */
/***************************** 6:2:91  OLD CODE **************
                        int origsize;
*******************************NEW CODE **********************/
    //int   origsize = 0;   /* sample type of file analysed */
    float F = 0.0f,
        F2 = 0.0f;      /* RWD: originally ~were~ float! */
    char    ch;     /* needed for crack (commandline interpreter)*/

    FILE    *fp;        /* auxiliary output file (-V option) */
/* new for PVOCEX */
    int insndfile = -1,outsndfile = -1;
    PSF_PROPS inprops,outprops;

    int pvfile = -1; /*in this version of pvoc, analysis file is input or output,
                     * but not both */
    int do_convert = 1;  /* set to 0 to write PVOC_COMPLEX instead */

    PVOCDATA *p_pvdata = NULL;
    WAVEFORMATEX *p_wfx = NULL;
    int flags,minheader = 0,do_messages = 1;
    float *analWindowBase;
    float *synWindowBase;
#ifdef USE_FFTW
    rfftwnd_plan forward_plan, inverse_plan;
    int in_fftw_size,out_fftw_size;
    float Ninv;
#endif

    float   beta = 6.8f,    /* parameter for Kaiser window */
        real,       /* real part of analysis data */
        imag,       /* imaginary part of analysis data */
        mag,        /* magnitude of analysis data */
        phase,      /* phase of analysis data */
        angleDif,   /* angle difference */
        RoverTwoPi, /* R/D divided by 2*Pi */
        TwoPioverR, /* 2*Pi divided by R/I */
        sum,        /* scale factor for renormalizing windows */
        ftot = 0.0f,    /* scale factor for calculating statistics */
        rIn,        /* decimated sampling rate */
        rOut,       /* pre-interpolated sampling rate */
        invR,       /* 1. / srate */
        time,       /* nI / srate */
        tvx0,       /* current x value of time-var function */
        tvx1,       /* next x value of time-var function */
        tvdx,       /* tvx1 - tvx0 */
        tvy0,       /* current y value of time-var function */
        tvy1,       /* next y value of time-var function */
        tvdy,       /* tvy1 - tvy0 */
        frac,       /* tvdy / tvdx */
        warp = 0.0f,    /* spectral envelope warp factor */
        R = 0.0f,       /* input sampling rate */
        P = 1.0f,       /* pitch scale factor */
        Pinv,       /* 1. / P */
        T = 1.0f;       /* time scale factor ( >1 to expand)*/

    int i,j,k,      /* index variables */
        Dd,     /* number of new inputs to read (Dd <= D) */
        Ii,     /* number of new outputs to write (Ii <= I) */
        N2,     /* N/2 */
        NO,     /* synthesis NO = N / P */
        NO2,        /* NO/2 */
        IO,     /* synthesis IO = I / P */
        IOi,        /* synthesis IOi = Ii / P */
        Mlen,
        Mf = 0,     /* flag for even M */
        Lf = 0,     /* flag for even L */
        Dfac,
        flag = 0,       /* end-of-input flag */
        C = 0,      /* flag for resynthesizing even or odd chans */
        Ci = 0,     /* flag for resynthesizing chans i to j only */
        Cj = 0,     /* flag for resynthesizing chans i to j only */
        CC = 0,     /* flag for selected channel resynthesis */
        V = 0,      /* verbose (summarize analysis) output flag */
        K = 0,      /* flag for Kaiser window */
        B = 0,      /* RWD : flag for rectancular (Box) window  */
        A = 0,      /* flag for analysis only */
        X = 0,      /* flag for magnitude output */
        E = 0,      /* flag for spectral envelope output */
        S = 0,      /* flag for synthesis only */
        tvflg = 0,  /* flag for time-varying time-scaling */
        tvnxt,      /* counter for stepping thru time-var func */
        tvlen;      /* length of time-varying function */       
    float   srate,      /* sample rate from header on stdin */
        arate;      /* sample rate for header on stdout if -A */
    float   timecheck = 0.0f;
    FUNCTION *tvpnt;    /* pointer to structure for time-var function*/
    int Nchans;         /* no of chans */
    /* RWD August 1 2006 */
    float* binfreqs = NULL;
    int vh = 0; /* do vonhann window */
/*RWD Oct12: write arrays to data files for examination! */
    
#ifdef NOTDEF
    int write_arrays = 0;
    int arraypos = 0;
    FILE* fpreal=NULL,*fpimag=NULL,*fpmag=NULL,*fpfreq=NULL;
    FILE* fpinmag = NULL, *fpinfreq = NULL;
    FILE* fpinPhase = NULL, *fpoutPhase = NULL;
    FILE* fpbinfreq =NULL, *fpbinphase = NULL;
    char fname[64] = {0};
    int write_fbin = 0;
    int write_pbin = 0;
    int thebin = 0;
    int binfreq = 1;
    int freq_out = 0;   // default: write input version of bin freq or phase
    int phase_out = 0;
    
#endif
    init_pvsys();
    psf_init();

    if(argc==2 && !strcmp(argv[1],"-h"))
        help();
    m = 0;
    for(n=1;n<argc;n++) {     
        if(m && *argv[n] !='-'){   /*don't count filenames starting with '-' !  */
            /*fprintf(stderr,"pvoc WARNING: filename starts with '-' ! \n"); */
            break;
        }
        if(*argv[n] =='-'){
            m++;
        }
    }
    if(argc-m < 3)                   /*RWD I don't trust this */
        usage();

/* call crack to interpret commandline */
    /*RWD NEW: B,K,m,v,x */
    /* and b,e,G,H */
    flags = 0;
    while ((ch = crack(argc,argv,"N|M|L|D|I|R|P|T|b|e|i|j|w|W|F|C|p|Q|H|XEABGSKVhmvx",0)) != 0){
        
        switch (ch) {
            case 'N': N    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'M': M    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'L': L    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'D': D    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'I': I    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'R': R    = (float) sfexpr(arg_option,1.0); break;
            case 'P': P    = (float) sfexpr(arg_option,1.0); break;
            case 'T': T    = (float) sfexpr(arg_option,1.0); break;
            case 'W': W    = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'F': F    = (float) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'C': C    = (int) (sfexpr(arg_option,1.0) + 0.5);CC=1; break;
            case 'i': Ci   = (int) (sfexpr(arg_option,1.0)+ 0.5); CC=1; break;
            case 'j': Cj   = (int) (sfexpr(arg_option,1.0)+ 0.5); CC=1; break;
            case 'b': nMin = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'e': nMax = (int) (sfexpr(arg_option,1.0) + 0.5); break;
            case 'w': 
                if (arg_option[0] == 0)
                    warp = -1.0f;
                else
                    warp = (float) sfexpr(arg_option,1.0);
                break;
/*RWD Oct12*/
#ifdef NOTDEF
            case 'p': 
                arraypos = (int) (sfexpr(arg_option,1.0) + 0.5); 
                write_arrays = 1; 
                break; 
            case 'H':
                write_pbin = 1;
                thebin = atoi(arg_option); /*(int) (sfexpr(arg_option,1.0) + 0.5);*/
                if(thebin < 0){
                    phase_out = 1;
                    thebin =  -thebin;                  
                }
                break;
                
            case 'Q':
                write_fbin = 1;
                thebin = atoi(arg_option); /*(int) (sfexpr(arg_option,1.0) + 0.5);*/
                if(thebin < 0){
                    freq_out = 1;
                    thebin = -thebin;                   
                }
                break;
#endif
            case 'A': 
                A = 1; break;
            case 'E':
                /*RWD*/
                fprintf(stderr,"\nSorry: -E option not yet supported for pvocex");
                exit(1);

                /*A = 1; E = 1; break;*/
            case 'X': 
                /*RWD*/
                fprintf(stderr,"\nSorry: -X option not yet supported for pvocex");
                exit(1);
                /*A = 1; X = 1; break;*/
            case 'S': 
                S = 1; break;
            case 'K':
                K = 1; 
                /*fprintf(stderr, "\npvoc: Kaiser window not implemented!\n");*/
                break;
            case 'G':
                vh = 1;
                break;
            case 'B':   /*RWD*/
                B = 1;
                break;
            case 'x':
                do_convert = 0;
                break;
            case 'V': 
                V = 1; 
                break;
            case 'm':           /*RWD new */
                minheader = 1;
                break;
            case 'v':
                do_messages = 0;
                break;
            case 'h': usage();  /* then exits normally */
            default: usage();   /* this exits with error */
        }
        flags++;
        //argv+= flags;
        //argc-= flags;
    }

    argv+= flags;
    argc-= flags;
    if(argc < 2){
        fprintf(stderr,"\ninsufficient arguments");
        usage();
        
    }

    /* RWD this message moved from below */
    if ((D != 0) && (T != 1.))
        fprintf(stderr,"pvocex: warning - don't specify both D and T\n");


    if(S){
        /* infile is analysis file*/
        p_pvdata = (PVOCDATA *) malloc(sizeof(PVOCDATA));
        if(p_pvdata==NULL){
            puts("\nNo Memory!");
            exit(1);
        }
        p_wfx = (WAVEFORMATEX *) malloc(sizeof(WAVEFORMATEX));
        if(p_wfx==NULL){
            puts("\nNo Memory!");
            exit(1);
        }


        pvfile  = pvoc_openfile(argv[1],p_pvdata,p_wfx);
        if(pvfile< 0){
            fprintf(stderr,"\npvocex: unable to open infile: %s",pvoc_errorstr());
            exit(1);
        }
        /* report format info */        
        printf("\nFormat info:");
        printf("\nChannels:        %d",p_wfx->nChannels);
        printf("\nSource Srate:    %d",p_wfx->nSamplesPerSec);
        printf("\nFrameCount:      %d",pvoc_framecount(pvfile));
        printf("\nFFT size:        %d",(p_pvdata->nAnalysisBins -1) *2);
        printf("\nWindow Size:     %d",p_pvdata->dwWinlen);
        printf("\nDecimation:      %d samples",p_pvdata->dwOverlap);
        printf("\nAnalysis Rate:   %.4f",p_pvdata->fAnalysisRate);
        printf("\nAnalysis Window: ");
        switch(p_pvdata->wWindowType){
        case(PVOC_HAMMING):
            printf("Hamming\n");
            break;
        case(PVOC_HANN):
            printf("Hanning\n");
            break;
        case(PVOC_KAISER):
            printf("Kaiser: param = %.4f\n",p_pvdata->fWindowParam);
            break;
        case(PVOC_CUSTOM):
            printf("Custom\n");
            break;
        default:
            printf("Default window\n");
            break;
        }
    }
    else{
        /* infile is soundfile */       
        insndfile = psf_sndOpen(argv[1],&inprops,0);
        if(insndfile < 0){
            fprintf(stderr,"\nUnable to open input soundfile %s",argv[1]);
            exit(1);
        }       
    }

    if (S) {
        /*wanted: 
         *  arate; Nchans -->N = (Nchans-1) * 2;
         *  Mlen --> M  (Winlen)
         *  Dfac --> D   (Overlap)
         *  srate-->R---srate = nSampleRate
         */
        arate = p_pvdata->fAnalysisRate;        
        Nchans = p_pvdata->nAnalysisBins;       
        if (N == 0)
            N = (Nchans - 1) * 2;
        Mlen = p_pvdata->dwWinlen;      
        if (M == 0)
            M = Mlen;
        Dfac = p_pvdata->dwOverlap;     
        if (D == 0)
            D = Dfac;
        srate = ((float) D * arate);
        if (R == 0.)
            R = srate;
        srate = R;
        
    }
    else {
        
        srate = (float) inprops.srate;  /* get input srate */
        if (R == 0.)
            R = srate;
        srate = R;
        
        Nchans = inprops.chans;
        /*RWD for pvocex this will no longer be a restriction, soon...*/
        if( Nchans != MONO){
            fprintf(stderr,"pvocex: input soundfile is not mono\n");
            exit(1);
        }
    }


#ifdef NOTDEF
    /*RWD: hidden until I can sort out the arglist mechanism */
    if (V) 
        if (argc > 2)
            fp = fopen(argv[3],"w");
    else
        fp = fopen("pvocex.s","w");

    if ((V != 1) && (arg_index < argc)) {
        tvpnt = gettsf(argv[arg_index++]);
        if(tvpnt==NULL) {       //RWD
            fprintf(stderr,"\npvocex: error opening function file %s\n",argv[arg_index-1]);
            exit(1);
        }
        tvflg = 1;
        tvx0 = (float)(tvpnt->fxval[0]);
        tvx1 = (float)(tvpnt->fxval[1]);
        tvy0 = (float)(tvpnt->fyval[0]);
        tvy1 = (float)(tvpnt->fyval[1]);

    fprintf(stderr,"name=%s\n", tvpnt->fname);
    fprintf(stderr,"type=%s\n", tvpnt->ftype);
    fprintf(stderr,"len=%d\n", tvpnt->flen);
    fprintf(stderr,"x0: %f  y0: %f  x1: %f  y1: %f\n",tvx0,tvy0,tvx1,tvy1);

        if (tvx0 != 0.){
    fprintf(stderr,"pvocex: warning - first x value in func must be 0\n");
            tvx0 = 0.0f;
        }
        if (tvy0 <= 0.){
    fprintf(stderr,"pvocex: invalid initial y value in time-vary func\n");
            exit(1);
        }
        tvdx = tvx1 - tvx0;
        if (tvdx <= 0.){
    fprintf(stderr,"pvocex: invalid x values in time-vary function\n");
            exit(1);
        }
        tvdy = tvy1 - tvy0;
        frac = tvdy / tvdx;
        tvnxt = 1;
        tvlen = tvpnt->flen;
    }
#endif
    /* RWD NB F now restored to float */
    if ((N != 0) && (F != 0.0f))
        fprintf(stderr,"pvocex: warning - don't specify both N and F\n");
    if ((N == 0) && (F == 0.0f))
        N = 256;
    else
        if (F != 0.0)
            N = (int)((float) R / F);
    N = N  + N%2;   /* Make N even */
    N2 = N / 2;

    F = /*(int)*/ ((float) R / N);
    F2 = F / 2.0f;

    if (W != -1) {
        if (M != 0)
            fprintf(stderr,"pvocex: warning - don't specify both M and W\n");
        else {
            switch(W){
                case 0: M = 4*N;
                    break;
                case 1: M = 2*N;
                    break;
                case 2: M = N;
                    break;
                case 3: M = N2;
                    break;
                default:
                    fprintf(stderr,"pvocex: warning - invalid W ignored\n");
                    break;
            }
        }
    }
    M = (M != 0 ? M : N );      /* RWD make double-window (W1) the default */
    Mf = 1 - M%2;

    L = (L != 0 ? L : M);
    Lf = 1 - L%2;

    if (M < 7)
        fprintf(stderr,"pvocex: warning - M is too small\n");

    ibuflen = 4 * M;
    obuflen = 4 * L;

    if (W == -1)
        W = (int)(3.322 * log10((double)(4. * N) / M));/* cosmetic */

    if ((A == 1) && (S == 1)) {
        fprintf(stderr,"pvocex: can't specify both -A and -S\n");
        exit(1);
        }

    if (Cj == 0) Cj = N2;
    if (Cj > N2) Cj = N2;
    if (Ci < 0) Ci = 0;

    // standard pitch-sclaing pvoc
    if(P != 1.0){
        if(T != 1.0)
            fprintf(stderr,"pvocex: warning - don't specify both T and P\n");

        if(tvflg){
            fprintf(stderr,"pvocex: can't specify -P with function\n");
            exit(1);
        }

        T = P;  /* pitch change is time change plus resamp */
    }
    //if ((D != 0) && (T != 1.))
    //  fprintf(stderr,"pvocex: warning - don't specify both D and T\n");

    if (T <= 0.){
        fprintf(stderr,"pvocex: invalid T = %f\n",T);
        exit(1);
    }

#ifdef OLD
    D = (D != 0 ? D : M/(8.0*max(1.0,T)) ); /* MCA fix - wrong types!! */
#endif
    D = (int)((D != 0 ? D : M/(8.0*(T > 1.0 ? T : 1.0))));

    if (D == 0){
        fprintf(stderr,"pvocex: warning - T greater than M/8 \n");
        D = 1;
    }

    I = (int)(I != 0 ? I : (float) T*D );

    if (I == 0){
        fprintf(stderr,"pvocex: T*D (or P*D) < 1 - increase M\n");
        exit(1);
    }

    if (tvflg){
        T = tvy0;   /* original T was maximum - used to set I */
        D = (int)((float) I / T);
        if (D < 1){
        fprintf(stderr,"pvocex: warning - can't expand by %f\n",T);
            D = 1;
        }
    }

    if (tvflg)
        if (warp != 0.)
            warp = T;   /* warp varies with T */

    T = ((float) I / D);
    if (P != 1.)
        P = T;

    NO = (int)((float) N / P);  /* synthesis transform will be NO points */
    NO = NO + NO%2;     /* make NO even */
    NO2 = NO / 2;
    P = ((float) N / NO);   /* ideally, N / NO = I / D = pitch change */
    Pinv = (float)(1.0/ P);

    if (warp == -1.)
        warp = P;
    if ((E == 1) && (warp == 0.))
        warp = 1.0f;


    if ((P != 1.) && (P != T))
         fprintf(stderr,"pvocex: warning P=%f not equal to T=%f\n",P,T);

    IO = (int)((float) I / P);
    nMax -= nMin;


    if (V) {
        fprintf(fp,"\nN: %d  M: %d  L: %d",N,M,L); 
        fprintf(fp,"  D: %d  I: %d  F: %f",D,I,F);
        fprintf(fp,"  R: %7.1f  P: %5.2f  T: %5.2f\n",R,P,T);
        if (CC)
            fprintf(fp,"C: %d    i: %d    j: %d\n",C,Ci,Cj);
        if (K)
            fprintf(fp,"---Kaiser Window---\n");
    }


/*RWD Oct12 generate text data files? */
#ifdef NOTDEF
    /* open data outfiles if requested */
    if(write_arrays){
        int err = 0;
        fpreal = fopen("inreal.dat","w+");
        if(fpreal==NULL) { err=1;goto fperr;}
        fpimag = fopen("inimag.dat","w+");
        if(fpimag==NULL) {err=2; goto fperr; }
        fpmag = fopen("outmag.dat","w+");
        if(fpmag==NULL) {err=3; goto fperr; }
        fpfreq = fopen("outfreq.dat","w+");
        if(fpfreq == NULL) {err=4; goto fperr;}
        fpinmag = fopen("inmag.dat","w+");
        if(fpinmag==NULL) {err=5; goto fperr; }
        fpinfreq = fopen("infreq.dat","w+");
        if(fpinfreq == NULL) {err=6; goto fperr;}
        fpinPhase = fopen("inphase.dat","w+");
        if(fpinPhase == NULL) {err=7; goto fperr;}
        fpoutPhase = fopen("outphase.dat","w+");
        if(fpoutPhase == NULL) {err=8; goto fperr;}
fperr:
        if(err){
            printf("error creating data file %d\n",err);
            return 1;
        }
    }
    if(write_fbin){         
        if(freq_out)    {
            sprintf(fname,"bin%doutfreq.dat",thebin);
            printf("writing output bin %d freq to %s.\n",thebin,fname);
        }
        else {
            sprintf(fname,"bin%dinfreq.dat",thebin);
            printf("writing input bin %d freq to %s.\n",thebin,fname);
        }
        fpbinfreq = fopen(fname,"w+");
        if(fpbinfreq==NULL){
            printf("error creating binfreq file %s\n",fname);
            return 1;
        }
    }

    if(write_pbin){
        if(phase_out) {
            sprintf(fname,"bin%doutphase.dat",thebin);
            printf("writing output bin %d phase to %s.\n",thebin,fname);
        }
        else {
            sprintf(fname,"bin%dinphase.dat",thebin);
            printf("writing input bin %d phase to %s.\n",thebin,fname);
        }
        fpbinphase = fopen(fname,"w+");
        if(fpbinphase==NULL){
            printf("error creating binphase file %s\n",fname);
            return 1;
        }
    }
#endif

    /*RWD: now, we can open outfile */
    if(A){
        /* create pvocex file - mono only, for now */
        pv_stype stype;
        pv_wtype wtype = PVOC_HAMMING;
        float winparam = 0.0f;
        if(K) {
            wtype = PVOC_KAISER;
            winparam = beta;
        }
        /* get this from infile - defaults to 16bit*/
        switch(inprops.samptype){
        
        case(PSF_SAMP_24):
            /* pvocex can in fact support both distinctly, using WAVE_EX, but I am being lazy*/
            stype = STYPE_24;
            break;
        case PSF_SAMP_32:
            stype = STYPE_32;
            break;
        case PSF_SAMP_IEEE_FLOAT:
            stype = STYPE_IEEE_FLOAT;
            break;
        default:
            stype = STYPE_16;
            break;
        }
        /* lots of default assumptions, for now */
        if(do_convert)
            pvfile  = pvoc_createfile(argv[2],N,D,1,PVOC_AMP_FREQ,(int)srate,stype,
                wtype,winparam,NULL,M);
        else
            pvfile  = pvoc_createfile(argv[2],N,D,1,PVOC_COMPLEX,(int)srate,stype,
                PVOC_RECT,winparam,NULL,M);
        if(pvfile < 0){
            fprintf(stderr,"\nUnable to create analysis file: %s",pvoc_errorstr());
            exit(1);
        }
    }
    else{
        /*  create soundfile: will have properties of infile
         *  i.e if infile is pvocex, the waveformatex fields of that determine
         *  the srate, sampletype, etc  of the outfile.
         */
        /*NB: if we have input soundfile, we can write outfile of the same type (WAVE, AIFF, etc)
         * but if we have pvoc input, need to decide outfile format from file extension if possible.
         * I am lazy, and just set WAVE format for now */
        int outchans =1,outsrate;
        psf_stype stype;
        

        if(pvfile >=0){
            outsrate = p_wfx->nSamplesPerSec;
            
            switch(p_wfx->wBitsPerSample){
            case(24):
                stype = PSF_SAMP_24;
                break;
            case(32):
                if(p_pvdata->wSourceFormat==WAVE_FORMAT_IEEE_FLOAT)
                    stype= PSF_SAMP_IEEE_FLOAT;
                else
                    stype = PSF_SAMP_32;
                break;          
            default:
                stype = PSF_SAMP_16;
                break;
            }

            /*check for window type - more to come...*/
            if(p_pvdata->wWindowType==PVOC_KAISER){
                K=1;
                if(p_pvdata->fWindowParam != 0.0f)
                    beta = p_pvdata->fWindowParam;
            }


        }
        else {
            
            outsrate = inprops.srate;
            stype = inprops.samptype;
        }

        if(insndfile >=0){
            outprops = inprops;
        /* keep all inprops  except... */
            outprops.srate = outsrate;
            outprops.samptype = stype;
        }
        else{
            outprops.chans = 1;
            outprops.srate = outsrate;
            outprops.samptype = stype;
            outprops.format = PSF_STDWAVE;
            outprops.chformat = STDWAVE;
        }

        
        /*last arg: no clipping of f/p samples: we use PEAK chunk, unless minheader */
        if((outsndfile = psf_sndCreate(argv[2],&outprops,0,minheader,PSF_CREATE_RDWR)) < 0){
            fprintf(stderr,"\nUnable to open outfile %s",argv[2]);
            exit(1);
        }
    }


#ifdef NOTDEF
    /* these formats not supported by pvocex (yet...)*/
    if (A){     /* Analysis only */
        if (E || X){
            Nchans = N2 + 1;            
            arate = ((float) Nchans * R / D);
        }
        else  {
            arate = ((float) R / D);
            Nchans = N + 2;         
        }       
    }
#endif  

    printf("\npvocex - analysis/synthesis beginning\n");
#ifdef NOTDEF
    if(thebin>=0){
        if(thebin > N2+1){
            printf("bin number %d too large.\n",thebin);
            return 1;
        }
    }
#endif
    /* set up analysis window: The window is assumed to be symmetric
        with M total points.  After the initial memory allocation,
        analWindow always points to the midpoint of the window
        (or one half sample to the right, if M is even); analWinLen
        is half the true window length (rounded down). Any low pass
        window will work; a Hamming window is generally fine,
        but a Kaiser is also available.  If the window duration is
        longer than the transform (M > N), then the window is
        multiplied by a sin(x)/x function to meet the condition:
        analWindow[Ni] = 0 for i != 0.  In either case, the
        window is renormalized so that the phase vocoder amplitude
        estimates are properly scaled.  The maximum allowable
        window duration is ibuflen/2. */


    analWindow = float_array(M+Mf);

    analWindowBase = analWindow;

    analWindow += (analWinLen = M/2);

    if (K){ 
        /*kaiser_(&M,analWindow,&analWinLen,(int*)1,&beta);*/
        /* RWD this function does whole window, not one half */
        kaiser(analWindowBase,/*analWinLen*/M,beta);            /*RWD */
        
    }
    else if(vh) {
         printf("Using Von Hann window\n");
         vonhann(analWindow,analWinLen,Mf);     
    }
    else if(B) {                                               /*RWD*/
        recwin(analWindow,analWinLen,Mf);
#ifdef _DEBUG
        printf("using rectangular window\n");
#endif
    }
    else
        hamming(analWindow,analWinLen,Mf);

    if(K==0)  /*RWD see above */
        for (i = 1; i <= analWinLen; i++)
            *(analWindow - i) = *(analWindow + i - Mf);

    if (M > N) {
        if (Mf)
        *analWindow *=(float)( (double)N * sin((double)PI*.5/N) /(double)( PI*.5));
        for (i = 1; i <= analWinLen; i++) 
            *(analWindow + i) *=(float)
            ((double)N * sin((double) (PI*(i+.5*Mf)/N)) / (PI*(i+.5*Mf)));   /* D.T. 2000*/
        for (i = 1; i <= analWinLen; i++)
            *(analWindow - i) = *(analWindow + i - Mf);
    }

    sum = 0.0f;
    for (i = -analWinLen; i <= analWinLen; i++)
        sum += *(analWindow + i);

    sum = (float)(2.0 / sum);       /*factor of 2 comes in later in trig identity*/

    for (i = -analWinLen; i <= analWinLen; i++)
        *(analWindow + i) *= sum;

    /* set up synthesis window:  For the minimal mean-square-error
        formulation (valid for N >= M), the synthesis window
        is identical to the analysis window (except for a
        scale factor), and both are even in length.  If N < M,
        then an interpolating synthesis window is used. */

    synWindow = float_array(L+Lf);
    synWindowBase = synWindow;
    synWindow += (synWinLen = L/2);
#ifdef USE_FFTW
    Ninv = (float) (1.0 / N);
#endif
    if (M <= N){
        if (K)  
            /*kaiser_(&M,synWindow,&synWinLen,(int*)1,&beta);*/
            kaiser(synWindowBase,/*synWinLen*/M,beta);          /*RWD */
        else if (B){                                    /*RWD*/
            recwin(synWindow,synWinLen,Lf);
#ifdef _DEBUG
            printf("using rectangular window\n");
#endif
        }
        else
            hamming(synWindow,synWinLen,Lf);
        if(K==0)
            for (i = 1; i <= synWinLen; i++)
                *(synWindow - i) = *(synWindow + i - Lf);

        for (i = -synWinLen; i <= synWinLen; i++)
            *(synWindow + i) *= sum;

        sum = 0.0f;
        for (i = -synWinLen; i <= synWinLen; i+=I)
            sum += *(synWindow + i) * *(synWindow + i);

        sum = (float)(1.0/ sum);
#ifdef USE_FFTW
        sum *= Ninv;
#endif
        for (i = -synWinLen; i <= synWinLen; i++)
            *(synWindow + i) *= sum;}
    else {  
        if (K)    /*RWD added: hope this is correct! */
            /*kaiser_(&M,synWindow,&synWinLen,(int*)1,&beta);*/
            kaiser(synWindowBase,/*synWinLen*/M,beta);          /*RWD */
        else

            hamming(synWindow,synWinLen,Lf);
        if(K==0)
            for (i = 1; i <= synWinLen; i++)
                *(synWindow - i) = *(synWindow + i - Lf);

        if (Lf)
            *synWindow *= (float)((double)IO * sin((double) (PI*.5/IO)) / (double)(PI*.5));
        for (i = 1; i <= synWinLen; i++) 
            *(synWindow + i) *=(float)
            ((double)IO * sin((double) (PI*(i+.5*Lf)/IO)) /(double) (PI*(i+.5*Lf)));
        for (i = 1; i <= synWinLen; i++)
            *(synWindow - i) = *(synWindow + i - Lf);

        sum = (float)(1.0/sum);
#ifdef USE_FFTW
        sum *= Ninv;
#endif
        for (i = -synWinLen; i <= synWinLen; i++)
            *(synWindow + i) *= sum;
    }

#ifdef USE_FFTW

    printf("Using FFTW\n");
    in_fftw_size = N;
    out_fftw_size = NO;
    forward_plan = rfftwnd_create_plan_specific(1,&in_fftw_size, 
        FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE,
        analWindowBase,1,NULL,1);
    inverse_plan = rfftwnd_create_plan_specific(1,&out_fftw_size, 
        FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE,
        synWindowBase,1,NULL,1);
    
#endif
    



    /* set up input buffer:  nextIn always points to the next empty
        word in the input buffer (i.e., the sample following
        sample number (n + analWinLen)).  If the buffer is full,
        then nextIn jumps back to the beginning, and the old
        values are written over. */

    input = float_array(ibuflen);

    nextIn = input;

    /* set up output buffer:  nextOut always points to the next word
        to be shifted out.  The shift is simulated by writing the
        value to the standard output and then setting that word
        of the buffer to zero.  When nextOut reaches the end of
        the buffer, it jumps back to the beginning.  */

    output =    float_array(obuflen);

    nextOut = output;
    /* set up analysis buffer for (N/2 + 1) channels: The input is real,
        so the other channels are redundant. oldInPhase is used
        in the conversion to remember the previous phase when
        calculating phase difference between successive samples. */

    anal =      float_array(N+2);
    banal = anal + 1;

    oldInPhase =    float_array(N2+1);
    maxAmp =    float_array(N2+1);
    avgAmp =    float_array(N2+1);
    avgFrq =    float_array(N2+1);
    env =       float_array(N2+1);
    /* set up synthesis buffer for (N/2 + 1) channels: (This is included
        only for clarity.)  oldOutPhase is used in the re-
        conversion to accumulate angle differences (actually angle
        difference per second). */

    syn =       float_array(NO+2);
    bsyn = syn + 1;

    oldOutPhase =   float_array(NO2+1);

    /* initialization: input time starts negative so that the rightmost
        edge of the analysis filter just catches the first non-zero
        input samples; output time is always T times input time. */

    
    if(R == 0) 
        exit(1);
    /* discard samples at start */
    /*for (i = 0; i < nMin; i++) */
        /*fgetfloat(nextIn,ifd);*/
#ifdef NOTDEF   
    if(psf_sndSeek(insndfile,nMin,PSF_SEEK_SET)){
        fprintf(stderr,"\nError seeking forward in infile");
        exit(1);
    }
#endif
    outCount = 0;
    rIn = ((float) R / D);
    rOut = ((float) R / I);
    invR =((float) 1. / R);
    RoverTwoPi = (float)(rIn / TWOPI);
    TwoPioverR = (float)(TWOPI / rOut);
    nI = -(analWinLen / D) * D; /* input time (in samples) */
    nO = (int)((float) T/P * nI);   /* output time (in samples) */
    Dd = analWinLen + nI + 1;   /* number of new inputs to read */
    Ii = 0;             /* number of new outputs to write */
    IOi = 0;
    flag = 1;

    stopwatch(1);

    /* main loop:  If nMax is not specified it is assumed to be very large
        and then readjusted when fgetfloat detects the end of input. */
    while(nI < (nMax + analWinLen)){
        int got;
        if(do_messages){
            time = nI * invR;
            if(time>timecheck){
                printf("\rInput time: %5.2f secs\t Output time: %5.2f secs",
                    time, nO * invR);
                fflush(stdout);
                timecheck = (float)(timecheck + TIME_INTERVAL);
            }
        }
        if (S == 1){            
            if((got = pvoc_getframes(pvfile,anal,1)) !=1){
                if(got < 0){
                    fprintf(stderr,"Error reading analysis file: %s",pvoc_errorstr());
                    exit(1);
                }
                else if(got == 0)
                    goto epilog;
            }
        }
        else {      /* prepare for analysis: read next Dd input values */
            {
                static float *sbuf = 0;
                static int sblen = 0;
                int got, tocp;
                float *sp;

                if(sblen < Dd) {
                    if(sbuf != 0)
                        free(sbuf);
                    if((sbuf = (float *)malloc(Dd*sizeof(float))) == 0) {
                        puts("pvocex: can't allocate float buffer\n");
                        exit(1);
                    }
                    sblen = Dd;
                }
                /*if((got = fgetsbuf(sbuf, Dd, ifd)) < 0) {*/
                if((got = psf_sndReadFloatFrames(insndfile,sbuf,Dd)) < 0){
                    fprintf(stderr,"\nError reading infile");
                    exit(1);
                }
                if(got < Dd)
                    Dd = got;
                sp = sbuf;

                tocp = min(got, input+ibuflen-nextIn);
                got -= tocp;
                while(tocp-- > 0)
                    *nextIn++ = *sp++;

                if(got > 0) {
                    nextIn -= ibuflen;
                    while(got-- > 0)
                        *nextIn++ = *sp++;
                }
                if (nextIn >= (input + ibuflen))
                    nextIn -= ibuflen;
            }

            if (nI > 0)
                for (i = Dd; i < D; i++){   /* zero fill at EOF */
                    *(nextIn++) = 0.0f;
                    if (nextIn >= (input + ibuflen))
                        nextIn -= ibuflen;
                }
    /* analysis: The analysis subroutine computes the complex output at
        time n of (N/2 + 1) of the phase vocoder channels.  It operates
        on input samples (n - analWinLen) thru (n + analWinLen) and
        expects to find these in input[(n +- analWinLen) mod ibuflen].
        It expects analWindow to point to the center of a
        symmetric window of length (2 * analWinLen +1).  It is the
        responsibility of the main program to ensure that these values
        are correct!  The results are returned in anal as succesive
        pairs of real and imaginary values for the lowest (N/2 + 1)
        channels.   The subroutines fft and reals together implement
        one efficient FFT call for a real input sequence.  */


        for (i = 0; i < N+2; i++) *(anal + i) = 0.0f;   /*initialize*/

        j = (nI - analWinLen - 1 + ibuflen) % ibuflen;  /*input pntr*/

        k = nI - analWinLen - 1;            /*time shift*/
        while (k < 0)
            k += N;
        k = k % N;
        for (i = -analWinLen; i <= analWinLen; i++) {
            if (++j >= ibuflen)
                j -= ibuflen;
            if (++k >= N)
                k -= N;
            *(anal + k) += *(analWindow + i) * *(input + j);
        }

#ifdef USE_FFTW
        rfftwnd_one_real_to_complex(forward_plan,anal,NULL);        
            //reals_(anal,banal,N2,-2);
#else
        fft_(anal,banal,1,N2,1,-2);
        reals_(anal,banal,N2,-2);
#endif
            
# ifdef NOTDEF
            /*RWD debug, test, learn, etc */
            /* print out real, imag arrays */
        if(write_arrays && outCount >= arraypos){
            /* fft real */
            int j;
            printf("\nwriting array info  from sample %d\n",outCount);
            assert(fpreal);
            assert(fpimag);
            
            for(j=0;j < N+2;j+=2){
                fprintf(fpreal,"%d\t%.6lf\n",j/2,anal[j]);
                fprintf(fpimag,"%d\t%.6lf\n",j/2,anal[j+1]);
            }
            /* after writing mag+freq; zet write_arrays to 0 so no duplication... */
        }
        
# endif


    /* conversion: The real and imaginary values in anal are converted to
        magnitude and angle-difference-per-second (assuming an 
        intermediate sampling rate of rIn) and are returned in
        anal. */

/*      i0 = anal - 2;
        i1 = anal - 1;
            oi = oldInPhase - 1;
        for (i = 0; i <= N2; i++)
*/

    if(do_convert){
        for(i=0,i0=anal,i1=anal+1,oi=oldInPhase; i <= N2; i++,i0+=2,i1+=2, oi++){
            real = *i0;
            imag = *i1;

            *i0 =(float) sqrt((double)(real * real + imag * imag));
        //  *i0 =  myhypot(real,imag);

                            /* phase unwrapping */
            if (*i0 == 0.)
                angleDif = 0.0f;

            else {

#ifdef USE_ATAN
                rratio = atan((double)imag/(double)real);               
                if(real<0.0) {
                    if(imag<0.0)
                        rratio -= PI;
                    else       
                        rratio += PI;
                }
/*
 * RWD: the modern form - sometimes faster, sometimes slower, depending ... 
 */
#else
                rratio = atan2((double)imag,(double)real);                              
#endif

                angleDif  = (phase = (float)rratio) - *oi;
                *oi = phase;
            }               
        /*  while */ if(angleDif >= PI)
                angleDif = (float)(angleDif - TWOPI);
        /*  while*/if (angleDif < -PI)
                angleDif = (float)(angleDif + TWOPI);

                        /* add in filter center freq.*/

            *i1 = angleDif * RoverTwoPi + ((float) i * F);

        }                    
    }  /* do_convert */

/* spectral envelope detection: this is a very crude peak picking algorithm
    which is used to detect and pre-warp the spectral envelope so that
    pitch transposition can be performed without altering timbre.
    The basic idea is to disallow large negative slopes between
    successive values of magnitude vs. frequency. */

        if (warp != 0.)
            warpse(anal,env,N2,warp);

        if (V) {
            ftot++;

            for (i = 0; i <= N2; i++){
                if (*(anal+2*i) > *(maxAmp+i))
                    *(maxAmp+i) = *(anal+2*i);
                *(avgAmp + i) += *(anal + 2*i);
                *(avgFrq + i) += *(anal + 2*i + 1);
            }
        }
    }

    /* if analysis only, write out interleaved instantaneous amplitudes
        and frequencies; otherwise perform resynthesis */

/**************  INSERTION_POINT FOR SPEC PROCESSES RECEIVING AND RETURNING
                    A SINGLE WINDOW 'anal' OF 'N2' CHANNELS *******************/

#ifdef NOTDEF
    /*RWD not supported yet: old CDP code here */
    if ((A==1) && (E==1))
        fputfbuf(env, N2+1, ofd);
        
    else if ((A == 1) && (X == 1)) {
        float *fp = anal;
        for (i=0; i <= N2; i++) {
            fputfloat(fp,ofd);
            fp += 2;
        }
    } 
    else  
#endif
    
        if (A == 1){
            /* sadly, a major bottle-neck on the PPC mac! */
            if(!pvoc_putframes(pvfile,anal,1)){
                fprintf(stderr,"\nError writing to analysis file: %s",pvoc_errorstr());
                exit(1);
            }
        }
        else {
                /* resynthesize only selected channels */
            if (CC){
                for (i = 0; i < Ci; i++)
                    *(anal+2*i) = 0.0f;
                for (i = Cj+1; i <= N2; i++)
                    *(anal+2*i) = 0.0f;
                if (C == 1)
                    for (i = Ci; i <= Cj; i++)
                        if (i%2 == 0)
                            *(anal+2*i) = 0.0f;
                if (C == 2)
                    for (i = Ci; i <= Cj; i++)
                        if (i%2 != 0)
                            *(anal+2*i) = 0.0f;
            }

        /* reconversion: The magnitude and angle-difference-per-second in syn
        (assuming an intermediate sampling rate of rOut) are
        converted to real and imaginary values and are returned in syn.
        This automatically incorporates the proper phase scaling for
        time modifications. */

            if (NO <= N){
                for (i = 0; i < NO+2; i++)
                    *(syn+i) = *(anal+i);
            }
            else {
                for (i = 0; i <= N+1; i++) 
                    *(syn+i) = *(anal+i);
                for (i = N+2; i < NO+2; i++) 
                    *(syn+i) = 0.0f;
            }

            for(i=0, i0=syn, i1=syn+1; i<= NO2; i++, i0+=2,  i1+=2){
                mag = *i0;
                //RWD NB no phase wrapping done here! 
                *(oldOutPhase + i) += *i1 - ((float) i * F);
                phase = *(oldOutPhase + i) * TwoPioverR;
                *i0 = (float)((double)mag * cos((double)phase));
                *i1 = (float)((double)mag * sin((double)phase));
            }
            if (P != 1.)
                for (i = 0; i < NO+2; i++)
                    *(syn+i) *= Pinv;

        /* synthesis: The synthesis subroutine uses the Weighted Overlap-Add
            technique to reconstruct the time-domain signal.  The (N/2 + 1)
            phase vocoder channel outputs at time n are inverse Fourier
            transformed, windowed, and added into the output array.  The
            subroutine thinks of output as a shift register in which 
            locations are referenced modulo obuflen.  Therefore, the main
            program must take care to zero each location which it "shifts"
            out (to standard output). The subroutines reals and fft
            together perform an efficient inverse FFT.  */

#ifdef USE_FFTW
            rfftwnd_one_complex_to_real(inverse_plan,(fftw_complex * )syn,NULL);

#else
            reals_(syn,bsyn,NO2,2);
            fft_(syn,bsyn,1,NO2,1,2);
#endif

            j = nO - synWinLen - 1;
            while (j < 0)
                j += obuflen;
            j = j % obuflen;

            k = nO - synWinLen - 1;
            while (k < 0)
                k += NO;
            k = k % NO;

            for (i = -synWinLen; i <= synWinLen; i++) { /*overlap-add*/
                if (++j >= obuflen)
                    j -= obuflen;
                if (++k >= NO)
                    k -= NO;

                *(output + j) += *(syn + k) * *(synWindow + i);

            }

            for (i = 0; i < IOi;){  /* shift out next IOi values */
                int j;
                int todo = min(IOi-i, output+obuflen-nextOut);              
                if(psf_sndWriteFloatFrames(outsndfile,nextOut,todo) < todo){
                    fprintf(stderr,"\nWarning: error writing data to outfile");
                    /* in case we can recover something anyway*/
                    goto epilog;
                }
                i += todo;
                outCount += todo;
                for(j = 0; j < todo; j++)
                    *nextOut++ = 0.0f;
                if (nextOut >= (output + obuflen))
                    nextOut -= obuflen;
            }
        }
                
        if (flag)   /* flag means do this operation only once */
            if ((nI > 0) && (Dd < D)){  /* EOF detected */
                flag = 0;
                nMax = nI + analWinLen - (D - Dd);
            }
/* time-varying time-scaling: linearly interpolate between x,y points */

        if (tvflg && (time > 0.)){
            while (tvflg && (time >= tvx1)) {
                if (++tvnxt >= tvlen)
                    tvflg = 0;
                else {
                    tvx0 = tvx1;
                    tvx1 = (float)(tvpnt->fxval[tvnxt]);
                    tvy0 = tvy1;
                    tvy1 = (float)(tvpnt->fyval[tvnxt]);
                    tvdx = tvx1 - tvx0;
                    if (tvdx <= 0.){
                        fprintf(stderr,"pvocex: invalid x values in time-vary function\n");
                        exit(2);
                    }
                    tvdy = tvy1 - tvy0;
                    frac = tvdy / tvdx;
                }
            }
            T = tvy0 + frac * (time - tvx0);
            if (T < ((float) 8. * I / M)){
                fprintf(stderr,"pvocex: warning - can't contract by %f\n",T);
                T = ((float) 8. * I / (M + 1));
            }
            D = (int)((double) I / T);
            if (D < 1){
                fprintf(stderr,"pvocex: warning - can't expand by %f\n",T);
                D = 1;
            }
            T = ((float) I / D);
            rIn = ((float) R / D);
            RoverTwoPi = (float)(rIn / TWOPI);
            if (warp != 0.)
                warp = T;
        }

    /*  D = some_function(nI);      for variable time-scaling */
    /*  rIn = ((float) R / D);      for variable time-scaling */
    /*  RoverTwoPi =  rIn / TwoPi;  for variable time-scaling */

        nI += D;                /* increment time */
        nO += IO;

    /* Dd = D except when the end of the sample stream intervenes */

        Dd = min(D, max(0, D+nMax-nI-analWinLen));


        if (nO > (synWinLen + I))
            Ii = I;
        else
            if (nO > synWinLen)
                Ii = nO - synWinLen;
            else {
                Ii = 0;
                for (i=nO+synWinLen; i<obuflen; i++)
                    if (i > 0)
                        *(output+i) = 0.0f;
            }
         IOi = (int)((float) Ii / P);
    }   /* End of main while loop */

    if (V) {
        if (ftot != 0.)
            ftot = (float)(1.0/ ftot);

        fprintf(fp,"\n   i         band       max amp     ");
        fprintf(fp,"avg amp    avg frq\n\n");
        for (i = 0; i <= N/2; i++) fprintf(fp,
            "%4d   %5.3f - %5.3f   %8.5f    %8.5f   %8.1f\n",
            2*i+1, i*F-F2, i*F+F2, *(maxAmp+i), *(avgAmp+i)*ftot, 
            *(avgFrq+i)*ftot);
        fprintf(fp,"\n");
    }

    if (A != 1) {
        nMaxOut = (int)((float) T / P * nMax);
        while (outCount <= nMaxOut){
            int todo = min(nMaxOut-outCount, output+obuflen-nextOut);
            if(todo == 0)
                break;          
            if(psf_sndWriteFloatFrames(outsndfile,nextOut,todo) < todo){
                fprintf(stderr,"Warning: error writing data to outfile\n");
                /* in case we can recover something anyway*/
                goto epilog;
            }
            outCount += todo;
            nextOut += todo;
            if (nextOut >= (output + obuflen))
                nextOut -= obuflen;
        }
    }
epilog:

    stopwatch(0);
    if(insndfile>=0) {
        psf_sndClose(insndfile);
        if(pvfile>=0){
            printf("Written %d analysis frames.\n",(int) pvoc_framecount(pvfile));
        }
    }
    if(outsndfile >=0)
        psf_sndClose(outsndfile);
    if(pvfile >=0)
        pvoc_closefile(pvfile);

    if(p_pvdata)
        free(p_pvdata);
    if(p_wfx)
        free(p_wfx);
    psf_finish();
#ifdef USE_FFTW
    rfftwnd_destroy_plan(forward_plan);
    rfftwnd_destroy_plan(inverse_plan);
#endif

/*RWD Oct12 */
#ifdef NOTDEF
    if(fpreal)
        fclose(fpreal);
    if(fpimag)
        fclose(fpimag);
    if(fpmag)
        fclose(fpmag);
    if(fpfreq)
        fclose(fpfreq);
    if(fpinmag)
        fclose(fpinmag);
    if(fpinfreq)
        fclose(fpinfreq);
    if(fpbinfreq)
        fclose(fpbinfreq);
    if(fpbinphase)
        fclose(fpbinphase);
#endif
#ifndef NOOVERCHK
    if(num_overflows > 0) {
        fprintf(stderr, "Warning: %d samples overflowed, and were clipped\n", num_overflows);
        fprintf(stderr, "\tmaximum sample was %f\n", maxsample);
        fprintf(stderr, "\tminimum sample was %f\n", minsample);
        if(-minsample > maxsample)
            maxsample = -minsample;
/*RWD: avoid div by zero...*/
        if(maxsample!=0.)
            fprintf(stderr, "You should reduce source level to avoid clipping: use gain of <=%lf\n", 1.0/maxsample);
    }
#endif
    return 0;
}


double
timer()
{
    struct timeb now;
    double secs, ticks; 
    ftime(&now);
    ticks = (double)now.millitm/(1000.0);
    secs = (double) now.time;

    return secs + ticks;
}


void            
stopwatch(int flag) 
{
    static double start;            
    int mins=0,hours=0;
    double end, secs=0;

    if(flag)
        start = timer();    
    else 
    {
        end    = timer(); 
        secs   = end - start;
        mins   = (int)(secs/60.0);
        secs  -= mins*60.0; 
        hours  = mins/60;
        mins  -= hours*60;

        fprintf(stderr,"\nElapsed time: ");
        if(hours > 0)
            fprintf(stderr,"%d h ", hours);
        if(mins > 0)
            fprintf(stderr,"%d m ", mins);
        fprintf(stderr,"%2.3lf s\n\n",secs);
    }
}
