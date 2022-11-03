/* ana2pvx.c: convert CDP analysis file to PVOC-EX */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sfsys.h>
/* need this to avoid clashes of GUID defs, etc */

#include "pvfileio.h"

#ifdef unix
int
_stricmp(const char *a, const char *b)
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
#else
// available in libsfsys
extern int _stricmp(const char *a, const char *b);
#endif

int main(int argc,char **argv)
{

    int i,anafd,pvxfd;
    int framesize, framesread = 0,numframes = 0;
    float *frame = NULL;
    char *ext;
    SFPROPS props;
    pv_stype stype;

    if(argc==1){
        printf("ANA2PVX: convert CDP analysis file to PVOC-EX file. \nv1.0 (c) RWD, CDP 2001.\n");
        printf("Usage:\tana2pvx inanalfile outfile.pvx\n");
        exit(1);
    }
    if(argc < 3){
        fprintf(stderr,"ana2pvx: insufficient arguments.\n"
                "Usage:\tana2pvx inanalfile outfile.pvx\n");
        exit(1);
    }

    ext = strrchr(argv[2],'.');
    
    if(ext==NULL || _stricmp(ext,".pvx")){
        fprintf(stderr,"anal2pvx: outfile name must use the .pvx extension.\n");
        exit(1);
    }

    sflinit("ana2pvx");
    init_pvsys();

    if((anafd = sndopenEx(argv[1],1,CDP_OPEN_RDONLY)) < 0){
        fprintf(stderr,"\nana2pvx: cannot open soundfile %s : %s",argv[1], rsferrstr);
        exit(1);
    }
    if(!snd_headread(anafd,&props)){
        fprintf(stderr,"\nana2pvx: error reading sfile header");
        fprintf(stderr,"\n%s",props_errstr);
        exit(1);
    }

    if(props.type != wt_analysis){
        fprintf(stderr,"anal2pvx: infile is not a CDP analysis file.\n");
        exit(1);
    }

    framesize = props.chans;
    frame = (float *) malloc(framesize * sizeof(float));
    if(frame==NULL){
        puts("ana2pvx: no memory!\n");
        exit(1);
    }
    switch(props.samptype){
    case FLOAT32:
        stype = STYPE_IEEE_FLOAT;
        break;
    case INT2424:
    case INT2432:
        stype = STYPE_24;
        break;
    case INT_32:
        stype = STYPE_32;
        break;
    default:
        stype = STYPE_16;
    }

    numframes = sndsizeEx(anafd);
    numframes /= framesize;
    if(numframes <= 0){
        fprintf(stderr,"ana2pvx: infile is empty!\n");
        exit(1);
    }

    pvxfd =  pvoc_createfile(argv[2], 
                 props.chans-2,props.decfac, 1,
                 PVOC_AMP_FREQ,props.origrate, 
                 stype,PVOC_HAMMING,0.0f,NULL,(DWORD) props.winlen);

    if(pvxfd < 0){
        fprintf(stderr,"ana2pvx: unable to open output file: %s\n",pvoc_errorstr());
        sndcloseEx(anafd);
        free(frame);
        exit(1);
    }

    for(i=0;i < numframes;i++){
        int got;
        if((got = fgetfbufEx(frame,framesize,anafd,0)) < framesize){
            if(got < 0){
                fprintf(stderr,"ana2pvx: error reading infile.\n");
                exit(1);
            }
            else if( got==0){
                fprintf(stderr,"ana2pvx: unexpected end of file after %d frames.\n",framesread);
                break;
            }
            else {
                fprintf(stderr,"ana2pvx:  incomplete frame read after %d frames: \n\tformat error in file?\n",framesread);
                break;
            }
        }
        if(!pvoc_putframes(pvxfd,frame,1)){
            fprintf(stderr,"ana2pvx: error writing outfile\n");
            exit(1);
        }
        framesread++;

        if(framesread%100==0)
            printf("copying: %d\r",framesread);
    }

    printf("ana2pvx: copied %d analysis frames to %s.\n",framesread, argv[2]);
    sndcloseEx(anafd);
    pvoc_closefile(pvxfd);
    free(frame);
    pvsys_release();
    return 0;
}
