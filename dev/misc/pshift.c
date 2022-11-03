#include <stdlib.h>
#include <assert.h>
#include "pshift.h"

#if defined _WIN32 && defined _MSC_VER
__inline static long round(double fval)
{
    int result;
    _asm{
        fld     fval
        fistp   result
        mov     eax,result
    }
    return result;
}

#else
static int round(double val)
{
    int k;
    k = (int)(fabs(val)+0.5);
    if(val < 0.0)
        k = -k;
    return k;
}
#endif



void get_amp_and_frq(const float *floatbuf,float *amp,float *freq,int clength)
{
    int cc, vc;
#ifdef _DEBUG
    assert(floatbuf != NULL);
    assert(amp != NULL);
    assert(freq != NULL);
#endif
    if(!(floatbuf && amp && freq))
        return;

    for( cc = 0 ,vc = 0; cc < clength; cc++, vc += 2){
        amp[cc]  = floatbuf[vc];
        freq[cc] = floatbuf[vc+1];
    }
    
}

void put_amp_and_frq(float *floatbuf,const float *amp, const float *freq,int clength)
{
    int cc, vc;
    if(!(floatbuf && amp && freq))
        return;

    for(cc = 0, vc = 0; cc < clength; cc++, vc += 2){
        floatbuf[vc] = amp[cc];
        floatbuf[vc+1] = freq[cc];
    }

}


void do_spectral_shiftp(float *amp, float *freq,double pitch,int clength)
{
    double shiftfac = pitch;
    int   j, k;
    
    /* RWD should never happen, but need to sort out Xtreme LFO mod to DC! */
    if(shiftfac==0.0)
        return;

    /* we do need to include the zero value somewhere...*/
    if( shiftfac >= 1.0) {
        j = clength-1;
        k  = round((double)j/shiftfac);
        while( k >= 0) {
            /*RWD*/
            if(j < 0)
                break;

            amp[j]  = amp[k];
            freq[j] = (float)(shiftfac * freq[k]);
            j-- ;
            k  = round((double)j/shiftfac);
        }
        for( k=j; k>= 0;k-- ) {      /*RWD was k > 0 */
            amp[k]  = 0.0f;
            freq[k] = 0.0f;
        }               
    } else if(shiftfac < 1.0){
        j = 0;
        k  = round((double)j/shiftfac);     
        while( k <= (clength-1)) {
            amp[j]  = amp[k];
            freq[j] = (float)(shiftfac * freq[k]);
            j++ ;
            k  = round((double)j/shiftfac);
        }
        for( k=j; k < clength; k++ ) {
            amp[k]  = 0.0f;
            freq[k] = 0.0f;
        }               
    }
}

