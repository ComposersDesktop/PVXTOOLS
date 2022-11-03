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
/* oscbank.h*/
#ifndef __OSCBANK_H_INCLUDED__
#define __OSCBANK_H_INCLUDED__


#ifndef MIN
#define MIN(x,y)    ( ((x)>(y)) ? (y) : (x) )
#endif

#ifndef MAX
#define MAX(x,y)    ( ((x)>(y)) ? (x) : (y) )
#endif

#define OSCBANK_TABLEN (16384)

class oscbank
{
public:
    oscbank();
    bool init(unsigned int maxoscs,float srate,unsigned int init_overlap);
    virtual ~oscbank();
    int synthesize_frame(float *pFrame,float *p_outbuf,float freqmod,unsigned int n_oscs);
    void set_overlap(unsigned int new_overlap) {
        overlap = new_overlap;
        one_over_overlap = 1.0f / overlap;
    }

protected:
    float tick(float *a,float *x,float *y){      /* must be inline */
        *x = *x - *a * *y;
        *y = *y + *a * *x;
        *y = MAX(-1.0f,*y);
        *y = MIN(1.0f,*y);
        return *y;
    }
    void clearall(void);
    void gensine(int tabsize);
private:
    
    float *a ,*x,*y;
    float *frame;  
    float *amps,*freqs,*indexes;
    float *lastamps,*lastfreqs;
    float *p_amp,*p_freq;
    float taboversr;
    float *accum;
    float *outbuf;
    float *p_outbuf;
    float maxamp,minamp;
    float minf_0,maxf_0;
    float pitchmod;
    float pi_over_sr;
    float one_over_overlap;
    float nyquist;
    unsigned int SampleRate,Channels;
    unsigned int framesize,buflen,overlap,winlen;
    unsigned int ob_size,tablen;
};


#endif
