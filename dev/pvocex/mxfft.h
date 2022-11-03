/*
 * Copyright (c) 1983-2013 Composers Desktop Project Ltd
 * http://www.composersdesktop.com
 * This file is part of the CDP System.
 * The CDP System is free software; you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The CDP System is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public
 * License along with the CDP System; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */
/* mxfft.h */
/* a proper header file for pvoc's famous fft routines! */



/*
        // RWD NB these are tiny arrays; where N=1024, maxf = 4, maxp = 6
        // where N=4096, maxf = 2, maxp = 8

    at = (float *) calloc(maxf,sizeof(float));
    ck = (float *) calloc(maxf,sizeof(float));
    bt = (float *) calloc(maxf,sizeof(float));
    sk = (float *) calloc(maxf,sizeof(float));
    np = (int *) calloc(maxp,sizeof(int));
*/

typedef struct mxfft_ 
{
    float* at;
    float* ck;
    float* bt;
    float* sk;
    int*    np;
    int ntot;

} MXFFTSTATE;

void fft2_(MXFFTSTATE* state, float* anal,float*  banal, int n, int N2, int nspn, int isn);
void reals_(float* anal,float* banal, int n,int isn);
MXFFTSTATE* create_mxstate(int nseg,int n, int nspn,int isn);
void free_mxstate(MXFFTSTATE* state);
