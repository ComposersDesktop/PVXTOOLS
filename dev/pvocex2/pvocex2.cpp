/*
 * Copyright (c) 2000,2022 Richard Dobson and Composers Desktop Project Ltd
 * http://www.rwdobson.com
 * http://www.composersdesktop.com
 *
 This file is part of the CDP System.
 
 The CDP System is free software; you can redistribute it
 and/or modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 The CDP System is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with the CDP System; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 02111-1307 USA
 *
 */
/* pvocex2.cpp*/
/*  demonstration stereo phase vocoder based on CARL PVOC
 *  This reads and writes all three defined frame formats,
 *  and supports a subset of the CARL flags.
 *  Now using portsf soundfile libary.
 *  
*/
/* RWD Feb 2008: fixed problem synthesising  analfile with v small number of frames */
 
#include <pvdefs.h>
extern "C"
{
#include <pvfileio.h>
#include <portsf.h>

}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>


#ifdef _DEBUG
#include <assert.h>
#endif

#ifdef unix
#include <ctype.h>
int stricmp(const char *a, const char *b);
int strnicmp(const char *a, const char *b, const int length);
#endif


#include "pvpp.h"
#include "oscbank.h"

#ifndef WAVE_FORMAT_IEEE_FLOAT
#define WAVE_FORMAT_IEEE_FLOAT (0x0003)
#endif

const int FFTLEN = 1024;



void stereo_split(float *inbuf,float *out_l,float *out_r,int insize);
void stereo_interl(float *in_l,float *in_r,float *out,int insize);
//int makeimpulsefile(const char *infilename,const char *outfilename);
//float *get_filterframe(const char *filename, int *pLength,int fnum);
//void complexmult(float *frame,const float *impulse,int length,int overlap);

#define DEFAULT_BUFLEN (32768)


void usage()
{
    printf("\n\nPVOCEX2 v0.15 RWD February 2008");
#ifdef USE_FFTW
    printf(": FFTW Version.");
#endif
    printf("\nStereo phase vocoder supporting fast convolution with impulse response.");
    printf("\nUsage: \n    pvocex2 [-A|-S][-F[fname]][-o[N]][-K|-H][-Px|-Tx][-Wx][-Nsamps][-m][-v]\n\t infile outfile"
           "\n-A      = perform Analysis only: output is .pvx file."
           "\n-S      = perform Synthesis only: input is .pvx file."
           "\n-o[Nosc]= use osc-bank resynthesis using Nosc oscillators (default = N/2 +1)"
           "\n-K      = use Kaiser window"
           "\n-H      = use von Hann window"
           "\n-R      = use Rectangular window"
           "\n          (default: Hamming window, or window from .pvx header)"
           "\n-Px     = Apply Pitch scaling rate x"
           "\n-Tx     = Apply Time scaling rate x"
           "\n          (NB: time/pitch scaling requires AMP_FREQ format)"    
           "\n-Wx     = set Window overlap factor (0,1,2,3, default 1)"
           "\n-Nsamps =  Set analysis window to <samps> samples (default 1024)"
           "\n-fx     = set Frame Type for analysis file."
           "\n          x = 0 = Amplitude,Frequency (default)"
           "\n          x = 1 = Amplitude,Phase ( = Soundhack format)"
           "\n          x = 2 = Complex (real,imaginary)"
           "\n         (NB: -f flag ignored unless -A is set.)"
           "\n-I      = show Format info of pvx file (no outfile required)."
           "\n-m      = write soundfile with minimum header (no PEAK data)"
           "\n-v      = suppress progress messages"
           "\nSoundfile formats supported: .wav.,.aiff.,.aif,.afc,.aifc\n");
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
        mins   = (long)(secs/60.0);
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


int main(int argc,char **argv)
{

    int fftlen = FFTLEN,windowsize  = FFTLEN,chans = 0,srate;
    int buflen = DEFAULT_BUFLEN,inbuflen,outbuflen,decfac,interp,inblocks;
    float nyquist, /*arate,*/ scalefac = 1.0f;
    
    phasevocoder *inpv_l = NULL,*outpv_l= NULL,*inpv_r=NULL,*outpv_r=NULL;
    pvoc_windowtype wtype = PVOC_HAMMING;
    pvoc_overlapfac ofac = PVOC_O_DEFAULT;
    pvoc_scaletype stype = PVOC_S_NONE;
    pvoc_frametype outframetype = PVOC_AMP_FREQ,inframetype = PVOC_AMP_FREQ;
    int insndfile = -1,outsndfile = -1;
    PSF_PROPS inprops,outprops;
    int i,j,T = 0,P=0,minheader = 0,do_timer = 1;
    int pvfile = -1; /* analysis file is input or output, but not both */
    
    PVOCDATA *p_pvdata = NULL;
    WAVEFORMATEX *p_wfx = NULL;
    oscbank *p_oscbank = NULL;
    oscbank *p_oscbank_r = NULL;
    unsigned int noscs = 0,ob_noscs_wanted = 0;
    float ob_freqfac = 1.0f /*,ob_timefac = 1.0f */;
    double oneovrsr;
    float *inbuf=NULL,*outbuf = NULL;
    float *inbuf_l = NULL,*outbuf_l = NULL;
    float *inbuf_r = NULL,*outbuf_r = NULL;
    float *pFrame_l = NULL,*pFrame_r = NULL;    
    bool do_Anal = true,do_Syn = true;  /* default is soundfile i/o */
    bool show_info = false;
    bool use_oscbank = false;

    char *ext;
    /* test params for spec pitch*/
    //float semitones = 11.5f;



    if(!init_pvsys()){
        puts("Unable to start pvsys.\n");
        return 1;
    }
    if(psf_init()){
        puts("unable to start portsf\n");
        return 1;
    }


    if(argc < 2){        /* possible -I infile.pvx */
        usage();
        return(1);
    }

    while(argc > 1 && argv[1][0]=='-'){
//        int userwinsize = 0;
        int wfac,fval;
        
        switch(argv[1][1]){
        case 'A':
            if(!do_Anal){
                fprintf(stderr,"\n-A cannot be used with -S.");
                return(1);
            }
            do_Syn = false;         
            break;
        case 'S':
            if(!do_Syn){
                fprintf(stderr,"\n-S cannot be used with -A.");
                return(1);
            }
            do_Anal = false;
            break;
        case 'N':
            if(argv[1][2]=='\0'){
                fprintf(stderr,"\nError: window size required for -N");
                usage();
                return(1);
            }
            fftlen = atoi(&(argv[1][2]));
            if(fftlen < 64){
                fprintf(stderr,"\nError: window size too small. Must be >= 64.");
                return(1);
            }
            fftlen = fftlen + fftlen%2;  /* make it even */
            break;
        case 'K':
            wtype = PVOC_KAISER;
            break;
        case 'H':
            wtype = PVOC_HANN;
            break;
        case 'R':
            wtype = PVOC_RECT;
            break;
        case 'W':
            if(argv[1][2]=='\0'){
                fprintf(stderr,"\n-W flag requires parameter.");
                usage();
                return(1);
            }
            wfac = atoi(&(argv[1][2]));
            switch(wfac){
            case 0:
                ofac = PVOC_O_W0;
                break;
            case 1:
                ofac = PVOC_O_W1;
                break;
            case 2:
                ofac = PVOC_O_W2;
                break;
            case 3:
                ofac = PVOC_O_W3;
                break;
            default:
                fprintf(stderr,"\n-W parameter out of range.");
                return(1);
            }
            break;
        case 'T':
            if(P){
                fprintf(stderr,"\n-T flag cannot be used with -P");
                return(1);
            }
            if(argv[1][2]=='\0'){
                fprintf(stderr,"\n-T flag requires parameter.");
                return(1);
            }
            scalefac = (float) atof(&(argv[1][2]));
            if(scalefac <= 0.0f){
                fprintf(stderr,"\nBad value for -T flag.");
                return(1);
            }
            stype = PVOC_S_TIME;
            T = 1;

            break;
        case 'P':
            if(T){
                fprintf(stderr,"\n-P flag cannot be used with -T");
                return(1);
            }
            if(argv[1][2]=='\0'){
                fprintf(stderr,"\n-P flag requires parameter.");
                return(1);
            }
            scalefac = (float) atof(&(argv[1][2]));
            if(scalefac <= 0.0f){
                fprintf(stderr,"\nBad value for -P flag.");
                return(1);
            }
            stype = PVOC_S_PITCH;
            P = 1;
            break;
        case 'f':
            if(argv[1][2]=='\0'){
                fprintf(stderr,"\n-F flag requires parameter.");
                return(1);
            }
            fval = atoi(&(argv[1][2]));
            if(fval < 0 || fval > 2){
                fprintf(stderr,"\n -F parameter out of range.");
                return(1);
            }
            outframetype = (pvoc_frametype) fval;
            break;
        case 'I':
            show_info = true;
            break;
        case 'm':
            minheader  =1;
            break;
        case 'v':
            do_timer = 0;
            break;
        case 'o':
            use_oscbank = true;
            if(argv[1][2]!='\0'){
                ob_noscs_wanted = atoi(&(argv[1][2]));
                if(ob_noscs_wanted <=0){
                    fprintf(stderr,"Error: for oscbank, at least 1 oscillator required.\n");
                    return(1);
                }
            }
            break;
        default:
            fprintf(stderr,"\nUnknown flag option %s ",argv[1]);
            return(1);
        }
        argc--;
        argv++;

    }
    /* 2 legal possibilities: infile and outfile, or -I used with infile only */
    if(argc< 2){
        fprintf(stderr,"\nInsufficient arguments.");
        usage();
        return(1);
    }
    if(argc < 3 && (argc < 2 && !show_info)){
        usage();
        return(1);
    }
    if(show_info && !do_Syn){
        fprintf(stderr,"\n-I flag cannot be used with -A");
        return(1);
    }
    if(show_info)
        do_Anal = false;
    /* warn about any ignored params, or bad extension: enforce .pvx here */
    if(!do_Syn){
        if(scalefac != 1.0) {
            if(P)
                printf("\nWarning: -P flag ignored for Analysis only.");
            else if(T)
                printf("\nWarning: -T flag ignored for analysis only.");
        }
        ext = strrchr(argv[2],'.');
        if(ext==NULL || stricmp(ext,".pvx")){
            fprintf(stderr,"\nError: outfile name must use extension .pvx");
            return(1);
        }
        if(use_oscbank){
            fprintf(stderr,"Warning: -o flag ignored for Analysis only.\n");
        }
    }



    if(do_Anal){
        /* infile is soundfile */       

        if((insndfile = psf_sndOpen(argv[1],&inprops, 0)) < 0){
            fprintf(stderr,"\nUnable to open input soundfile %s",argv[1]);
            return(1);
        }

        chans = inprops.chans;
        /*RWD for pvocex this will no longer be a restriction, soon...*/
        if( chans > 2){
            fprintf(stderr,"\nSorry, can only read mono or stereo soundfiles at present.\n");
            return(1);
        }
        srate = inprops.srate;
        if(srate <=0){
            fprintf(stderr,"\nbad srate found: corrupted file?\n");
            return(1);
        }
        oneovrsr = 1.0 / (double) srate;
        inpv_l = new phasevocoder();
        if(inpv_l==NULL){
            puts("\nunable to create inpv_l");
            return(1);
        }

        /* use generic init function */
        if(!inpv_l->init(srate,fftlen,ofac,scalefac,stype,wtype,PVPP_OFFLINE)){
            fprintf(stderr,"\nUnable to initialize inpv_l.");
            return(1);
        }

        if(chans==2){
            inpv_r = new phasevocoder();
            if(inpv_r==NULL){
                puts("\nunable to create inpv_r");
                return(1);
            }       
            if(!inpv_r->init(srate,fftlen,ofac,scalefac,stype,wtype,PVPP_OFFLINE)){
                fprintf(stderr,"\nUnable to initialize inpv_r.");
                return(1);
            }
        }
        decfac = inpv_l->anal_overlap();
        windowsize = inpv_l->winlen();
    }
    else{
        /*input is pvx file */
        char *ext = strrchr(argv[1],'.');
        if(ext==NULL || stricmp(ext,".pvx")){
            fprintf(stderr,"\nInfile does not have .pvx extension.");
            return(1);
        }

        p_pvdata = (PVOCDATA *) malloc(sizeof(PVOCDATA));
        if(p_pvdata==NULL){
            puts("\nNo Memory!");
            return(1);
        }
        p_wfx = (WAVEFORMATEX *) malloc(sizeof(WAVEFORMATEX));
        if(p_wfx==NULL){
            puts("\nNo Memory!");
            return(1);
        }

        pvfile  = pvoc_openfile(argv[1],p_pvdata,p_wfx);
        if(pvfile< 0){
            fprintf(stderr,"\npvocex2: unable to open infile: %s",pvoc_errorstr());
            return(1);
        }

        /*stick to mono/stereo for now */
        if(p_wfx->nChannels > 2){
            fprintf(stderr,"\nSorry: can only read mono or stereo pvx files.");
            return(1);
        }
        chans = p_wfx->nChannels;
        fftlen = (p_pvdata->nAnalysisBins - 1) *2;
        windowsize = p_pvdata->dwWinlen;
        decfac = p_pvdata->dwOverlap;
        inframetype = (pvoc_frametype) p_pvdata->wAnalFormat;
        srate = p_wfx->nSamplesPerSec;
        if(srate <=0){
            fprintf(stderr,"\nbad srate found: corrupted file?\n");
            return(1);
        }
        oneovrsr = 1.0 / (double) srate;

        wtype = (pvoc_windowtype) p_pvdata->wWindowType;

        /* report format info */
        if(show_info){
            printf("\nFormat info:");
            printf("\nChannels:          %d",chans);
            printf("\nSource Srate:      %d",p_wfx->nSamplesPerSec);
            printf("\nSource sampletype: ");
            printf("%d bit ",p_wfx->wBitsPerSample);
            if(p_wfx->wBitsPerSample ==32)
                if(p_pvdata->wSourceFormat==WAVE_FORMAT_IEEE_FLOAT)
                    printf("floats");
            printf("\nFrameCount:        %d",pvoc_framecount(pvfile));
            printf("\nFFT size:          %d",fftlen);
            printf("\nWindow Size:       %d",windowsize);
            printf("\nDecimation:        %d samples",decfac);
            printf("\nAnalysis Rate:     %.4f",p_pvdata->fAnalysisRate);
            printf("\nAnalysis Window:   ");
        
            switch(wtype){
            case(PVOC_HAMMING):
                printf("Hamming");
                break;
            case(PVOC_HANN):
                printf("von Hann");
                break;
            case(PVOC_KAISER):
                printf("Kaiser: param = %.4f",p_pvdata->fWindowParam);
                break;
            case(PVOC_CUSTOM):
                printf("Custom");
                break;
            case(PVOC_RECT):
                printf("Rectangular");
                break;
            default:
                printf("Default window");
                break;
            }
            printf("\nAnalysis Type:     ");
            switch(inframetype){
            case(PVOC_AMP_FREQ):
                printf("Amplitude-Frequency\n");
                break;
            case(PVOC_AMP_PHASE):
                printf("Amplitude-Phase\n");
                break;
            case(PVOC_COMPLEX):
                printf("Complex (Real-Imaginary)\n");
                break;
            }
            if(argc < 3)
                return(0);
        }
        
        /* time/pitch scaling only supported for amp-freq at the moment */
        if(scalefac != 1.0 && inframetype!= PVOC_AMP_FREQ)
            printf("\nWarning: time/pitch rescaling only supported for AMP_FREQ format.");
        /* create phasevocoders from .pvx properties */
        inpv_l = new phasevocoder();
        if(inpv_l==NULL){
            puts("\nunable to create inpv_l");
            return(1);
        }
        /* use specific init function */
        
        if(!inpv_l->init(srate,fftlen,windowsize,decfac,scalefac,stype,wtype,PVPP_OFFLINE)){
            fprintf(stderr,"\nUnable to initialize inpv_l.");
            return(1);
        }

        if(chans==2){
            inpv_r = new phasevocoder();
            if(inpv_r==NULL){
                puts("\nunable to create inpv_r");
                return(1);
            }       
            if(!inpv_r->init(srate,fftlen,windowsize,decfac,scalefac,stype,wtype,PVPP_OFFLINE)){
                fprintf(stderr,"\nUnable to initialize inpv_r.");
                return(1);
            }
        }
    }

    if(do_Syn){
        /* input is pvx file, or from analysis, output is sfile */
        outpv_l = new phasevocoder();
        if(outpv_l==NULL){
            puts("\nunable to create outpv");
            delete inpv_l;
            return(1);
        }
        if(p_pvdata){
            if(!outpv_l->init(srate,fftlen,windowsize,decfac,scalefac,stype,wtype,PVPP_OFFLINE)){
                fprintf(stderr,"\nUnable to initialize outpv_l.");
                delete outpv_l;
                return(1);
            }

        }       
        else if(!outpv_l->init(srate,fftlen,ofac,scalefac,stype,wtype,PVPP_OFFLINE)){
            fprintf(stderr,"\nUnable to initialize outpv_l.");
            delete outpv_l;
            return(1);
        }
        /* or another possibility is to init directly from inpv_l */

        if(chans==2){       
            outpv_r = new phasevocoder();
            if(outpv_r==NULL){
                puts("\nunable to create outpv");
                delete inpv_r;
                return(1);
            }
            if(p_pvdata){
                if(!outpv_r->init(srate,fftlen,windowsize,decfac,scalefac,stype,wtype,PVPP_OFFLINE)){
                    fprintf(stderr,"\nUnable to initialize outpv_l.");
                    delete outpv_l;
                    return(1);
                }

            }       
            else if(!outpv_r->init(srate,fftlen,ofac,scalefac,stype,wtype,PVPP_OFFLINE)){
                fprintf(stderr,"\nUnable to initialize outpv_r.");
                delete outpv_r;
                return(1);
            }
        }

        /*create output sfile */
        psf_stype stype;
        psf_format outformat;
        /* will it be aiff, wav, etc? */
        outformat = psf_getFormatExt(argv[2]);
        if(outformat==PSF_FMT_UNKNOWN){
            fprintf(stderr,"\nOutfile name has unrecognized extension.");
            return(1);
        }
        if(pvfile >=0){
            /* get sfile sampletype from pvx file */
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
#ifdef _DEBUG
                stype = PSF_SAMP_IEEE_FLOAT;
#else
                stype = PSF_SAMP_16;
#endif
                break;
            }
        }
        else            
            stype = inprops.samptype;

        /*the one gotcha: if output is floats and format is aiff, change to aifc */
        if(stype==PSF_SAMP_IEEE_FLOAT && outformat==PSF_AIFF){
            fprintf(stderr,"Warning: AIFF output written as AIFC for float samples\n");
            outformat = PSF_AIFC;
        }

        outprops = inprops;
        outprops.chans  = chans;
        outprops.srate = srate;
        outprops.format = outformat;
        outprops.samptype = stype;
        outprops.chformat = STDWAVE;
        /*last arg: no clipping of f/p samples: we use PEAK chunk */
        if((outsndfile = psf_sndCreate(argv[2],&outprops,0,minheader,PSF_CREATE_RDWR)) <0 ){
            fprintf(stderr,"\nUnable to open outfile %s",argv[2]);
            return(1);
        }
    }
    else{
        /* output is .pvx file, input is sfile */
        
#ifdef _DEBUG
        assert(insndfile);
#endif
        pv_stype stype;
                
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
        pvfile  = pvoc_createfile(argv[2],fftlen,decfac,chans,outframetype,srate,stype,
            wtype,0.0,NULL,windowsize);
        if(pvfile < 0){
            fprintf(stderr,"\nUnable to create analysis file: %s",pvoc_errorstr());
            return(1);
        }
    }

    decfac = inpv_l->anal_overlap();
    if(do_Syn)
        interp = outpv_l->syn_interp();
    else
        interp = decfac;
    /* try to set buffer lengths to suit decfac,interp blocksizes.
     * we can get smaller increments than these, but never larger, so these sizes are safe.*/

    /*make sure buflen is at least fftlen */
    if(fftlen > buflen)
        buflen = fftlen;

    inbuflen = buflen + (buflen % decfac);
    inblocks = inbuflen / decfac;
    inbuflen = inblocks * decfac;

    outbuflen = inblocks  * interp;
    
    if(p_pvdata && show_info)       /* add to Format report */
        printf("\nInterpolation:     %d",interp);
#ifdef NOTDEF
    printf("\nUsing ");
    switch(wtype){
    case(PVOC_HAMMING):
        printf("Hamming Window");
        break;
    case(PVOC_HANN):
        printf("von Hann Window");
        break;
    case(PVOC_KAISER):
        printf("Kaiser Window");
        break;
    case(PVOC_CUSTOM):
        printf("Custom Window");
        break;
    default:
        printf("Default Window");
        break;
    }
    printf("\nAnalysis Window Size: %d",windowsize);
#endif  
    /* main i/o buffers */
    if(do_Anal){
        inbuf = new float[inbuflen * chans];
        inbuf_l = new float[inbuflen];
        if(chans==2)
            inbuf_r = new float[inbuflen];
    }
    if(do_Syn){     
        outbuf = new float[outbuflen * chans];
        outbuf_l = new float[outbuflen];        
        if(chans==2){
            outbuf_r = new float[outbuflen];            
        }
    }

    /* channel buffers */
    
    pFrame_l = new float[fftlen+2];
    if(chans==2)                
        pFrame_r = new float[fftlen+2];
    
    printf("\n");
    unsigned int nbins = (fftlen + 2) / 2;   /* i.e clength  */
    nyquist = (float) srate * 0.5f;
    //float  chwidth = nyquist/(float)(nbins - 1);


    if(do_Syn && use_oscbank) {
        if(ob_noscs_wanted==0)
            noscs = nbins;
        else
            noscs = min(ob_noscs_wanted,nbins);
        p_oscbank = new oscbank();
        if(!p_oscbank->init(noscs,(float) srate,interp)){
            fprintf(stderr,"Unable to create oscbank.\n");
            return(1);
        }
        if(chans==2){
            p_oscbank_r = new oscbank();
            if(!p_oscbank_r->init(noscs,(float) srate,interp)){
                fprintf(stderr,"Unable to create oscbank.\n");
                return(1);
            }
        }
        if(P)
            ob_freqfac = scalefac;
        printf("Using Oscillator-bank resynthesis, with %d oscillators.\n",noscs);
    }

    stopwatch(1);

    long written,thisblock,sampsread;
    long samps_to_write= 0;
    double intime= 0.0,outtime = 0.0;
    if(do_Anal && do_Syn){
        if(chans==1){
            while((sampsread  = psf_sndReadFloatFrames(insndfile,inbuf,inbuflen))  > 0){                
                written = thisblock = samps_to_write = 0;
                /* zeropad to full buflen */
                if(sampsread < inbuflen){
                    for(i = sampsread;i< inbuflen;i++)
                        inbuf[i] = 0.0f;
                    sampsread = inbuflen;
                }
                for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                    inpv_l->generate_frame(inbuf+i,pFrame_l,decfac,(pvoc_frametype)PVOC_AMP_FREQ);

                    /*****  INSERT FX HERE! ******/

                    if(use_oscbank)
                        thisblock = p_oscbank->synthesize_frame(pFrame_l,outbuf+j,ob_freqfac,noscs);
                    else
                        thisblock = outpv_l->process_frame(pFrame_l,outbuf+j,(pvoc_frametype)PVOC_AMP_FREQ);
                    samps_to_write += thisblock;
                }
#ifdef _DEBUG
                if(samps_to_write > outbuflen)
                    assert(0);
#endif
                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,samps_to_write)) != samps_to_write){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    intime += (double)sampsread * oneovrsr;
                    outtime += (double)written * oneovrsr;              
                    printf("Input time: %.2lf\t Output time: %.2lf\r",intime,outtime);
                }
            }

            /* write out remainder */
            written = thisblock = samps_to_write = 0;
            sampsread = fftlen;
            for(i = 0;i < sampsread;i++)
                inbuf[i] = 0;

            for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                    inpv_l->generate_frame(inbuf+i,pFrame_l,decfac,(pvoc_frametype)PVOC_AMP_FREQ);

                    /*****  INSERT FX HERE! ******/
                    if(use_oscbank)
                        thisblock = p_oscbank->synthesize_frame(pFrame_l,outbuf+j,ob_freqfac,noscs);
                    else
                        thisblock = outpv_l->process_frame(pFrame_l,outbuf+j,(pvoc_frametype)PVOC_AMP_FREQ);
                    samps_to_write += thisblock;
                }
#ifdef _DEBUG
                if(samps_to_write > outbuflen)
                    assert(0);
#endif
                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,samps_to_write)) != samps_to_write){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    intime += (double)sampsread * oneovrsr;
                    outtime += (double)written * oneovrsr;              
                    printf("Input time: %.2lf\t Output time: %.2lf\r",intime,outtime);
                }
        }
        else {      /* Anal+Synth: stereo */
            while((sampsread = psf_sndReadFloatFrames(insndfile,inbuf,inbuflen)) > 0){              
                written = thisblock = samps_to_write = 0;
                /* zeropad to full buflen */
                if(sampsread < inbuflen){
                    for(i = sampsread*chans;i< inbuflen*chans;i++)
                        inbuf[i] = 0.0f;
                    sampsread = inbuflen;
                }
                stereo_split(inbuf,inbuf_l,inbuf_r,sampsread*chans);
                
                for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                    inpv_l->generate_frame(inbuf_l+i,pFrame_l,decfac,(pvoc_frametype)PVOC_AMP_FREQ);
                    inpv_r->generate_frame(inbuf_r+i,pFrame_r,decfac,(pvoc_frametype)PVOC_AMP_FREQ);
                    /*****  INSERT TRANSFORM HERE! ******/ 
                    if(use_oscbank){
                        thisblock = p_oscbank->synthesize_frame(pFrame_l,outbuf_l+j,ob_freqfac,noscs);
                        thisblock = p_oscbank_r->synthesize_frame(pFrame_r,outbuf_r+j,ob_freqfac,noscs);

                    }
                    else{
                        thisblock = outpv_l->process_frame(pFrame_l,outbuf_l+j,(pvoc_frametype)PVOC_AMP_FREQ);
                        thisblock = outpv_r->process_frame(pFrame_r,outbuf_r+j,(pvoc_frametype)PVOC_AMP_FREQ);
                    }
                    samps_to_write += thisblock;
                }           
#ifdef _DEBUG
                if(samps_to_write * chans > outbuflen * chans)
                    assert(0);
#endif
                stereo_interl(outbuf_l,outbuf_r,outbuf,samps_to_write);

                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,samps_to_write)) != samps_to_write){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    intime += (double)sampsread * oneovrsr;
                    outtime += (double)written * oneovrsr; 
                    printf("Input time: %.2lf\t Output time: %.2lf\r",intime,outtime);
                }
            }

            /* write out remainder */
            written = thisblock = samps_to_write = 0;
            sampsread = fftlen;
            for(i=0;i < sampsread * chans;i++)
                inbuf[i] = 0;
            stereo_split(inbuf,inbuf_l,inbuf_r,sampsread*chans);
                
            for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                inpv_l->generate_frame(inbuf_l+i,pFrame_l,decfac,(pvoc_frametype)PVOC_AMP_FREQ);
                inpv_r->generate_frame(inbuf_r+i,pFrame_r,decfac,(pvoc_frametype)PVOC_AMP_FREQ);
                /*****  INSERT TRANSFORM HERE! ******/ 
                if(use_oscbank){
                    thisblock = p_oscbank->synthesize_frame(pFrame_l,outbuf_l+j,ob_freqfac,noscs);
                    thisblock = p_oscbank_r->synthesize_frame(pFrame_r,outbuf_r+j,ob_freqfac,noscs);
                }
                else{
                    thisblock = outpv_l->process_frame(pFrame_l,outbuf_l+j,(pvoc_frametype)PVOC_AMP_FREQ);
                    thisblock = outpv_r->process_frame(pFrame_r,outbuf_r+j,(pvoc_frametype)PVOC_AMP_FREQ);
                }
                samps_to_write += thisblock;
            }           
#ifdef _DEBUG
            if(samps_to_write * chans > outbuflen * chans)
                assert(0);
#endif
            stereo_interl(outbuf_l,outbuf_r,outbuf,samps_to_write);

            if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,samps_to_write)) != samps_to_write){
                fprintf(stderr,"\nerror writing to outfile");
                return(1);
            }
            if(do_timer){
                intime += (double)sampsread * oneovrsr;
                outtime += (double)written * oneovrsr; 
                printf("Input time: %.2lf\t Output time: %.2lf\r",intime,outtime);
            }
        }
    }
    else if(do_Anal){
#ifdef _DEBUG
        assert(pvfile >= 0);
        assert(insndfile);
#endif
        if(chans==1){           
            while((sampsread  = psf_sndReadFloatFrames(insndfile,inbuf,inbuflen))  > 0){
                written = thisblock = samps_to_write = 0;
                /* zeropad to full buflen - actually just for final source buffer */
                if(sampsread < inbuflen){
                    for(i = sampsread;i< inbuflen;i++)
                        inbuf[i] = 0.0f;
                    sampsread = inbuflen;
                }
                for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                    inpv_l->generate_frame(inbuf+i,pFrame_l,decfac,outframetype);
                    if(!pvoc_putframes(pvfile,pFrame_l,1)){
                        fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                        return(1);
                    }
                }
                if(do_timer){
                    intime += (double)sampsread * oneovrsr;
                    printf("Input time: %.2lf\r",intime);
                }
            }
            /* write out remaining frames */
            sampsread = fftlen;
            for(i = 0;i< sampsread;i++)
                inbuf[i] = 0.0f;
            for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                inpv_l->generate_frame(inbuf+i,pFrame_l,decfac,outframetype);
                if(!pvoc_putframes(pvfile,pFrame_l,1)){
                    fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                    return(1);
                }
            }
            if(do_timer){
                intime += (double)sampsread * oneovrsr;
                printf("Input time: %.2lf\r",intime);
            }
        }
        else {         // doAnal: stereo
            while((sampsread = psf_sndReadFloatFrames(insndfile,inbuf,inbuflen)) > 0){
                written = thisblock = samps_to_write = 0;
                /* zeropad to full buflen */
                if(sampsread < inbuflen){
                    for(i = sampsread*chans;i< inbuflen*chans;i++)
                        inbuf[i] = 0.0f;
                    sampsread = inbuflen;
                }
                stereo_split(inbuf,inbuf_l,inbuf_r,sampsread*chans);
                
                for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                    inpv_l->generate_frame(inbuf_l+i,pFrame_l,decfac,outframetype);
                    inpv_r->generate_frame(inbuf_r+i,pFrame_r,decfac,outframetype);
                    if(!pvoc_putframes(pvfile,pFrame_l,1)){
                        fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                        return(1);
                    }
                    if(!pvoc_putframes(pvfile,pFrame_r,1)){
                        fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                        return(1);
                    }
                }
                if(do_timer){
                    intime += (double)sampsread * oneovrsr;
                    printf("Input time: %.2lf\r",intime);
                }
            }

            /* write out remaining frames */
            sampsread = fftlen * chans;
            for(i = 0;i< sampsread;i++)
                inbuf[i] = 0.0f;
            for(i=0,j=0;i < sampsread; i+= decfac,j+= thisblock){
                inpv_l->generate_frame(inbuf_l+i,pFrame_l,decfac,outframetype);
                inpv_r->generate_frame(inbuf_r+i,pFrame_r,decfac,outframetype);
                if(!pvoc_putframes(pvfile,pFrame_l,1)){
                    fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                    return(1);
                }
                if(!pvoc_putframes(pvfile,pFrame_r,1)){
                    fprintf(stderr,"\nError writing analysis frames: %s",pvoc_errorstr());
                    return(1);
                }
            }
            if(do_timer){
                intime += (double)sampsread * oneovrsr;
                printf("Input time: %.2lf\r",intime);
            }
        }
    }
    else {     /* must be Synth only! */
#ifdef _DEBUG
        assert(do_Syn);
        assert(outsndfile>=0);
        assert(pvfile >= 0);
        assert(p_pvdata);
#endif
        int numframes = pvoc_framecount(pvfile);
        int got = 0,framesread = 0;
        int j =0 ;
        int blocks = 0;   /* variable retval from process_frame, so must count blocks
                                to ensure we don't overrun buffers */
        if(chans==1){           
            while(pvoc_getframes(pvfile,pFrame_l,1)){
                if(p_oscbank)
                    got = p_oscbank->synthesize_frame(pFrame_l,outbuf+j,ob_freqfac,noscs);
                else
                    got = outpv_l->process_frame(pFrame_l,outbuf+j,inframetype);
                framesread++;
                j+= got;
                if(++blocks < inblocks)
                    continue;
#ifdef _DEBUG
                if(j > outbuflen)
                    assert(0);
#endif
                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    outtime += (double)written * oneovrsr;
                    printf("Output time: %.2lf\r",outtime);
                }
                /* reset buffer counters for next outbufs-worth*/
                j = 0;
                blocks = 0;
                
            }
            /* write out anything left in this block, before doing pvoc remainder*/
            if(j > 0){
                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    outtime += (double)written * oneovrsr;
                    printf("Output time: %.2lf\r",outtime);
                }
            }
            
            /* generate remaining blocks  */
            j =0;
            while(j < fftlen){
                //got = outpv_l->process_frame(pFrame_l,outbuf+j,inframetype);
                if(p_oscbank)                               
                    got = p_oscbank->synthesize_frame(pFrame_l,outbuf+j,ob_freqfac,noscs);
                else
                    got = outpv_l->process_frame(pFrame_l,outbuf+j,inframetype);
                j+= got;
#ifdef _DEBUG
                if(j > outbuflen)
                    assert(0);
#endif
            }
            if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                fprintf(stderr,"\nerror writing to outfile");
                return(1);
            }
            if(do_timer){
                outtime += (double)written * oneovrsr;
                printf("Output time: %.2lf\r",outtime);
            }
        }
        else {          
            while(pvoc_getframes(pvfile,pFrame_l,1)){
                 /* no second frame: must be error! */
                if(!pvoc_getframes(pvfile,pFrame_r,1)){
                    fprintf(stderr,"\nError reading frame for Ch 2: %s", pvoc_errorstr());
                    return(1);
                }
                if(use_oscbank){
                    got = p_oscbank->synthesize_frame(pFrame_l,outbuf_l+j,ob_freqfac,noscs);
                    got = p_oscbank_r->synthesize_frame(pFrame_r,outbuf_r+j,ob_freqfac,noscs);
                }
                else{
                    got = outpv_l->process_frame(pFrame_l,outbuf_l+j,inframetype);
                    got = outpv_r->process_frame(pFrame_r,outbuf_r+j,inframetype);
                }
                framesread+=2;
                j += got;
                blocks++;
                if(blocks < inblocks)
                    continue;

#ifdef _DEBUG
                if(j  > outbuflen)
                    assert(0);
#endif
                stereo_interl(outbuf_l,outbuf_r,outbuf,j);

                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    outtime += (double)written * oneovrsr;
                    printf("Output time: %.2lf\r",outtime);
                }
                j = 0;
                blocks = 0;
            }
            if(j> 0){
                stereo_interl(outbuf_l,outbuf_r,outbuf,j);
                if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                    fprintf(stderr,"\nerror writing to outfile");
                    return(1);
                }
                if(do_timer){
                    outtime += (double)written * oneovrsr;
                    printf("Output time: %.2lf\r",outtime);
                }
            }

            /* generate remaining blocks */

            j =0;
            while(j < fftlen){
                if(use_oscbank){
                    got = p_oscbank->synthesize_frame(pFrame_l,outbuf_l+j,ob_freqfac,noscs);
                    got = p_oscbank_r->synthesize_frame(pFrame_r,outbuf_r+j,ob_freqfac,noscs);
                }
                else{
                    got = outpv_l->process_frame(pFrame_l,outbuf+j,inframetype);
                    got = outpv_r->process_frame(pFrame_r,outbuf_r+j,inframetype);
                }
                j+= got;
#ifdef _DEBUG
                if(j > outbuflen)
                    assert(0);
#endif
            }
            stereo_interl(outbuf_l,outbuf_r,outbuf,j);
            if((written = psf_sndWriteFloatFrames(outsndfile,outbuf,j)) != j){
                fprintf(stderr,"\nerror writing to outfile");
                return(1);
            }
            if(do_timer){
                outtime += (double)written * oneovrsr;
                printf("Output time: %.2lf\r",outtime);
            }
        }
        if(framesread < numframes)
            fprintf(stderr,"\nWarning: read only %d frames out of %d.",framesread,numframes);
    }
    

    stopwatch(0);
    if(insndfile >=0)
        psf_sndClose(insndfile);
    if(outsndfile >=0)
        psf_sndClose(outsndfile);
    
    if(pvfile >= 0)
        pvoc_closefile(pvfile);

    pvsys_release();
    psf_finish();

    if(p_oscbank)
        delete p_oscbank;
    if(p_oscbank_r)
        delete p_oscbank_r;
    
    if(inbuf)
        delete [] inbuf;
    if(outbuf)
        delete [] outbuf;
    if(inbuf_l)
        delete [] inbuf_l;
    if(outbuf_l)
        delete [] outbuf_l;
    if(inpv_l)
        delete inpv_l;
    if(outpv_l)
        delete outpv_l;
    if(pFrame_l)
        delete [] pFrame_l;
    if(chans==2){
        if(inbuf_r)
            delete [] inbuf_r;
        if(outbuf_r)
            delete [] outbuf_r;
        if(inpv_r)
            delete inpv_r;
        if(outpv_r)
            delete outpv_r;
        if(pFrame_r)
            delete [] pFrame_r; 
    }
    return 0;
}

void stereo_split(float *inbuf,float *out_l,float *out_r,int insize)
{
    long i;
    float *pfl_l,*pfl_r,*pfl_i;

    pfl_i = inbuf;
    pfl_l = out_l;
    pfl_r = out_r;
    /* still need to confirm this is actually faster than an increasing index! */
    for(i = insize/2;i;--i){
        *pfl_l++ = *pfl_i++;
        *pfl_r++ = *pfl_i++;
    }

}

void stereo_interl(float *in_l,float *in_r,float *out,int insize)
{
    long i;
    float *pfl_l,*pfl_r,*pfl_o;
    pfl_o = out;
    pfl_l = in_l;
    pfl_r = in_r;

    for(i=insize;i;--i){
        *pfl_o++ = *pfl_l++;
        *pfl_o++ = *pfl_r++;
    }
}

void complexmult(float *frame,const float *impulse,int length, int overlap)
{
    float re,im;
    float scalefac;
    int i,j;
    /* need some scaling, but is this calc reasonable? */
    scalefac  = (float) (1.0 / pow(2,(double) length / (double)overlap));
    

    for(i=0,j = 1;i < length;i+=2,j+=2){
        re = frame[i] * impulse[i] - frame[j] * impulse[j];
        im = frame[i] * impulse[j] + frame[j]* impulse[i];
        frame[i] = re * scalefac;
        frame[j] = im * scalefac;
    }
}

#ifdef unix
int stricmp(const char *a, const char *b)
{
    while(*a != '\0' && *b != '\0') {
        int ca = islower(*a) ? toupper(*a) : *a;
        int cb = islower(*b) ? toupper(*b) : *b;
        
        if(ca < cb)
            return -1;
        if(ca > cb)
            return 1;
        
        a++;
        b++;
    }
    if(*a == '\0' && *b == '\0')
        return 0;
    if(*a != '\0')
        return 1;
    return -1;
}

int
strnicmp(const char *a, const char *b, const int length)
{
    int len = length;
    
    while(*a != '\0' && *b != '\0') {
        int ca = islower(*a) ? toupper(*a) : *a;
        int cb = islower(*b) ? toupper(*b) : *b;
        
        if(len-- < 1)
            return 0;
        
        if(ca < cb)
            return -1;
        if(ca > cb)
            return 1;
        
        a++;
        b++;
    }
    if(*a == '\0' && *b == '\0')
        return 0;
    if(*a != '\0')
        return 1;
    return -1;
}
#endif
