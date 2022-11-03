/* pvconv.c
 * convert Csound or Soundhack pvoc file to  PVOC_EX file 
 * Initial Version 0.1 RWD 25:5:2000: all rights reserved: work in progress! 
 * many assumptions, chiefly that the input format is correct!
 */


#include <stdio.h>
#include <stdlib.h>


#include <pvfileio.h>

#ifndef PI
#define PI	(3.14159265358979323846f)
#endif
#ifndef TWOPI
#define TWOPI (2.0f * PI)
#endif
#define PVMAGIC 517730	/* look at it upside-down, esp on a 7-seg display */
#define ERBEMAGIC 0x45726265				 /*"Erbe" */
#define PVDFLTBYTS 4
#define PVAMPFAC (1.0f / 256.0f)

#define REVDWBYTES(t)	( (((t)&0xff) << 24) | (((t)&0xff00) << 8) | (((t)&0xff0000) >> 8) | (((t)>>24) & 0xff) )

typedef struct pvstruct
    {
    long	magic;			/* magic number to identify */
    long	headBsize;		/* byte offset from start to data */
    long	dataBsize;		/* number of bytes of data */
    long	dataFormat;	       	/* (int) format specifier */
    float	samplingRate;		/* of original sample */
    long	channels;		/* (int) mono/stereo etc */
    long 	frameSize;		/* size of FFT frames (2^n) */
    long	frameIncr;		/* # new samples each frame */
    long	frameBsize;		/* bytes in each file frame */
    long	frameFormat;		/* (int) how words are org'd in frms */
    float	minFreq;		/* freq in Hz of lowest bin (exists) */
    float	maxFreq;		/* freq in Hz of highest (or next) */
    long	freqFormat;		/* (int) flag for log/lin frq */
    char	info[PVDFLTBYTS];	/* extendable byte area */
    } PVSTRUCT;

enum pv_filetype {PV_CSOUND,PV_SOUNDHACK};

void convert_to_ampfreq(float *frame,float *phase,long size, float fund, float Pifac,float nyquist);

void main (int argc, char **argv)
{

	int ifd,ofd,chans;
	int i;
	int numframes, framesize;           
	int need_bytereverse = 0;
	long framesread = 0;
	long overlap,winlen;
	float fundamental, *pframe = NULL,*phase= NULL;
	float pvampfac = 1.0f,RovrTwoPi;
	double maxamp_frame,maxamp,arate;
		
    FILE *fp ;    
    PVSTRUCT *phdr;
    
	enum pv_filetype pvtype = PV_CSOUND;	/* assume */
	
	PVOCDATA pvdata;
	WAVEFORMATEX fmtex;




    phdr =  (PVSTRUCT *)malloc( sizeof(PVSTRUCT)) ;
    
	printf("\nPVCONV: convert Csound or Soundhack pvoc file to .pvx");
	printf("\n        Initial Version 0.1 RWD 25:5:2000");	   

	if(argc < 3){
		fprintf(stderr,"\n\nUsage: pvconv [-r] infile.pv outfile.pvx");
		fprintf(stderr,"\n   -r: byte-reverse infile data\n");
		exit(1);
	}

	while(argv[1][0]=='-'){
		switch(argv[1][1]){
		case 'r':
			need_bytereverse = 1;
			break;
		default:
			fprintf(stderr,"\nunknown flag option %s",argv[1]);
			exit(1);
		}
		argc--;
		argv++;
	}
	if(argc < 3){
		fprintf(stderr,"Usage: pvconv [-r] infile.pv outfile.pvx\n");
		exit(1);
	}

	if ( ( fp = fopen( argv[1], "rb") ) ==NULL ) {
      fprintf( stderr, "pvconv: Unable to open '%s'\n Does it exist?",
               argv[1] ) ;
      exit( -1 ) ;
    }

	if(fread(phdr, 1, sizeof(PVSTRUCT), fp) != sizeof(PVSTRUCT)){
		fprintf(stderr,"\nerror reading pv file header.");
		exit(1);
	}

	if(need_bytereverse){
		int i,n = sizeof(PVSTRUCT) / sizeof(long);
		long *pl = (long *) phdr;
		for(i=0; i < n;i++)
			pl[i] = REVDWBYTES(pl[i]);

	}

	switch(phdr->magic){
	case PVMAGIC:			  /* assume it is type 7 */
		
		break;
	case ERBEMAGIC:
		printf("\nreading Soundhack file");
		if(phdr->frameFormat==3)
			pvtype = PV_SOUNDHACK;
		else if(phdr->frameFormat==7)		
			pvtype = PV_CSOUND;
		break;
	default:    
		fprintf( stderr, "'%s' is not a recognized pvoc file\n", argv[1] ) ;
		exit( -1 ) ;
		break;
    }

    framesize = phdr->frameSize+2;
    numframes =((phdr->dataBsize/sizeof(float)) / phdr->frameSize);
	chans = phdr->channels;
	pframe = (float *) malloc(framesize* sizeof(float));
	if(pframe==NULL){
		puts("\nNo memory for frame buffer!");
		exit(1);
	}
	for(i=0;i < framesize;i++)
		pframe[i] = 0.0f;

	if(pvtype==PV_SOUNDHACK){
		phase =  (float *) malloc((framesize >> 1) * sizeof(float));
		if(phase==NULL){
			puts("\nNo memory for Soundhack phase buffer!");
			exit(1);
		}
		for(i=0;i < (framesize>>1);i++)
			phase[i] = 0.0f;

	}
	

	
	fundamental = phdr->samplingRate / (float) (phdr->frameSize);
		
	winlen = phdr->frameSize;
	overlap = phdr->frameIncr;
	arate = (float)(phdr->samplingRate / (float)overlap);

	RovrTwoPi = (float)(arate / TWOPI);
	maxamp = 0.0;

	init_pvsys();

	/* alternatively, can write Soundhack file directly as PVOC_AMP_PHASE */
	/* assume hamming window: no info in Csound or Soundhack file */
	ofd = pvoc_createfile(argv[2],phdr->frameSize,overlap,chans,PVOC_AMP_FREQ,(long)(phdr->samplingRate),
		STYPE_16,PVOC_HAMMING,0.0f,NULL,(DWORD) winlen);
	if(ofd < 0){
		fprintf(stderr,"\nunable to open outfile: %s",pvoc_errorstr());
		exit(1);
	}
		
	if(pvtype==PV_CSOUND)
		pvampfac = 1.0f / (float) phdr->frameIncr;	
	else
		pvampfac = 1.0;	


	printf("\nreading %d-channel file....\n",chans);
	if(pvtype==PV_SOUNDHACK){
		/* write empty first frame */
		if(!pvoc_putframes(ofd,pframe,1)){
		    fprintf(stderr,"\nError writing frame to outfile");
			exit(1);
		}
	}
	if(need_bytereverse){
		while(fread(pframe,1,framesize* sizeof(float),fp) == framesize* sizeof(float)){
		
			long *pl = (long *) pframe;
			maxamp_frame = 0.0;
			for(i=0;i < framesize;i +=2) {				
			
				pl[i] = REVDWBYTES(pl[i]);
				pl[i+1] = REVDWBYTES(pl[i+1]);
				pframe[i] *= pvampfac;
				maxamp_frame = max(maxamp_frame,(double)(pframe[i]));
			}
			if(pvtype == PV_SOUNDHACK && phdr->frameFormat==3){
				convert_to_ampfreq(pframe,phase,framesize, fundamental, RovrTwoPi,phdr->samplingRate/2.0f);
				
			}
#ifdef _DEBUG
			printf("\nFrame %d: maxamp = %.6lf",framesread,maxamp_frame);
#endif

			if(!pvoc_putframes(ofd,pframe,1)){
			    fprintf(stderr,"\n%s: %s",argv[2],pvoc_errorstr());
				exit(1);
			}
			framesread++;
			maxamp = max(maxamp,maxamp_frame);
		}



	}

	else{
		while(fread(pframe,1,framesize* sizeof(float),fp) == framesize* sizeof(float)){
			maxamp_frame = 0.0;
			for(i=0;i < framesize;i +=2) {				
				pframe[i] *= pvampfac;
				maxamp_frame = max(maxamp_frame,(double)(pframe[i]));
			}
			if(pvtype == PV_SOUNDHACK)
				convert_to_ampfreq(pframe,phase,framesize, fundamental, RovrTwoPi,phdr->samplingRate/2.0f);
#ifdef _DEBUG
				printf("\nFrame %d: maxamp = %.6lf",framesread,maxamp_frame);
#endif
			if(!pvoc_putframes(ofd,pframe,1)){
			    fprintf(stderr,"\nError writing frame to outfile");
				exit(1);
			}

			framesread++;
		}
	}
	
	if(!pvoc_closefile(ofd)){
		fprintf(stderr,"\nError closing PVOC-EX file");
		exit(1);
	}

	if(framesread != numframes)
		printf("\nWARNING: infile has %d frames, written %d frames",numframes, framesread);
	printf("\nConverted %d frames, maxamp = %.4lf",framesread,maxamp);
	fclose(fp);
	if(pframe)
		free(pframe);
	if(phase)
		free(phase);
	free(phdr);

	/*  TEST: read file in, and check header */

	ifd = pvoc_openfile(argv[2],&pvdata,&fmtex);
	if(ifd < 0){
		fprintf(stderr,"\nUnable to reopen file: %s",pvoc_errorstr());
	}
	else{
		printf("\nfile reopened to read header:");
		printf("\nFormat info:");
		printf("\nChannels: %d",fmtex.nChannels);
		printf("\nSource Samplerate: %d",fmtex.nSamplesPerSec);
		printf("\nFrameCount: %d",pvoc_framecount(ifd));
		printf("\nWindow Size: %d",pvdata.dwWinlen);
		printf("\nDecimation: %d samples",pvdata.dwOverlap);
		printf("\nAnalysis Rate: %.4f",pvdata.fAnalysisRate);
	}
	pvoc_closefile(ifd);
	pvsys_release();

}

/* this will probably go in pvfileio.c eventually */

void convert_to_ampfreq(float *frame,float *phase, long size, float fund, float Pifac,float Nyquist)
{
	int i;
	long nbins;
	float frq,angleDif,thisphase;
	for(i=0,nbins = 0; i < size; i += 2,nbins++){

		thisphase = frame[i+1];		
			angleDif = thisphase - phase[nbins];

		phase[nbins ] = thisphase;
		while (angleDif > PI)
			angleDif -= TWOPI;
		while (angleDif < -PI)
			angleDif += TWOPI;


		frq = angleDif * Pifac + ((float)nbins * fund);
		
		frame[i+1]  = frq;
#ifdef DEBUG_VERBOSE
				
		printf("\n[%d]:\tAMP = %.6f\tPH = %.6f\tANG = %.6f\tFRQ = %.6f",
			i/2,frame[i],thisphase,angleDif,frq);
		
#endif
		


	}
}
