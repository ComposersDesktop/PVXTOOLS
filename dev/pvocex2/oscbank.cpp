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
/* oscbank.cpp*/
#include <math.h>
#include <memory.h>
#ifdef _DEBUG
#include <assert.h>
#endif
#include "oscbank.h"

#ifndef PI
#define PI  (3.14159265358979323846)
#endif

#define TWOPI   (2*PI)
#ifndef NULL
#define NULL (0)
#endif


oscbank::oscbank()
{
    a = x = y = NULL;
    frame = NULL;  
    amps = freqs = indexes =NULL;
    lastamps = lastfreqs = NULL;
    accum = NULL;
    outbuf = NULL;
    p_outbuf = NULL;    
    pitchmod = 1.0; 
    ob_size = 0;
}

oscbank::~oscbank()
{
    clearall();

}




bool oscbank::init(unsigned int maxoscs,float srate,unsigned int init_overlap)
{
    double oneovrsr;
    try {
    
        a           = new float[maxoscs];
        x           = new float[maxoscs];
        y           = new float[maxoscs];
        amps        = new float[maxoscs];
        lastamps    = new float[maxoscs];
        freqs       = new float[maxoscs];
        lastfreqs   = new float[maxoscs];
        indexes     = new float[maxoscs];

    }
    catch(...){
        clearall();

        return false;
    }

    for(unsigned int i=0;i < maxoscs;i++){
        a[i] = 0.0f;
        x[i] = 1.0f;
        y[i] = 0.0f;
        amps[i] = lastamps[i] = freqs[i] = indexes[i] = 0.0f;
    }


    taboversr = (float)tablen / srate;
    overlap = init_overlap;
    one_over_overlap = 1.0f / overlap;
    oneovrsr = 1.0 / srate;
    pi_over_sr = (float)(PI * oneovrsr);
    nyquist = srate * 0.5f;
    ob_size = maxoscs;
    return true;
}


int oscbank::synthesize_frame(float *pFrame,float *p_outbuf,float freqmod, 
                               unsigned int n_oscs)
{
    unsigned int i,j;

#ifdef _DEBUG
    assert(ob_size);
    assert(n_oscs <= ob_size);
    assert(pFrame);
    assert(p_outbuf);
#endif
    if(n_oscs ==0 || ob_size==0)
        return 0;
    /* clear output buffer */
    memset(p_outbuf,0,overlap * sizeof(float));

    /*update amps, freqs*/

    for(i=0;i < n_oscs;i++){
        amps[i] =   pFrame[i*2];
        freqs[i] = (float)(fabs)(pFrame[(i*2)+1]) * freqmod;
        if(freqs[i] >= nyquist)
            a[i] = 0.0f;
        else
            a[i] = (float )(2.0 * sin(freqs[i] * pi_over_sr));
    }
    /* we can interp amplitude, but freq is not so easy with this method...*/
    for(i=0;i < n_oscs;i++) {
        float thisamp = lastamps[i];
        float delta_amp = (amps[i] - thisamp) * one_over_overlap;
        
        for(j=0;j < overlap;j++) {
            p_outbuf[j] += thisamp *  tick(a+i,x+i,y+i);
            thisamp += delta_amp;
        }
        lastamps[i] = amps[i];
    }

    return overlap;
}


void oscbank::clearall(void)
{
    if(a) {
        delete [] a;
        a = 0;
    }
    if(x) {
        delete [] x;
        x = 0;
    }
    if(y) {
        delete [] y;
        y = 0;
    }
    if(amps) {
        delete [] amps;
        amps = 0;
    }
    if(lastamps) {
        delete [] lastamps;
        lastamps = 0;
    }
    if(freqs) {
        delete [] freqs;
        freqs = 0;
    }
    if(lastfreqs) {
        delete [] lastfreqs;
        lastfreqs = 0;
    }
    if(indexes) {
        delete [] indexes;
        indexes = 0;
    }

}
