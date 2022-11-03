/* cmdutils.h*/
/*
 * Copyright (c) 1984,2022, Composers Desktop Project Ltd
 
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
/* command-line parsing: support for POLISH, CRACK, etc. etc */
#define LRPN 200            /* longest possible expression */
#define UNOPS "{-}"         /* unary operators for polish() */ 
#define BINOPS "{^,%}{*,/}{+,-}"    /* binary operators for polish() */
#define POSTOPS "{dB,K,k,S,s,sec,secs,m,ms,L,C,R}" /* postops for polish() */
#define EXUNOPS "{sin,cos,atan,ln,exp,floor,abs,rand,sqrt}{-}"
#define EXPOSTOPS "{dB,K,k,Deg,invs,MM}"
#define SFPOSTOPS "{dB,K,k,S,s,m,ms}"




typedef struct func {
    char   *ftype;    /* function type: MONO_IN_X, or XY_PAIRS */
    char   *fname;    /* function name */
    long    flen;     /* function length */
    double *fxval;    /* x function values */
    double *fyval;    /* y function values */
} FUNCTION;


char  *polish(char *expression,char *unops,char *binops,char *postops);
double sfexpr(char *string,double timefac);
int    sffield(char **input,char *string,char *iglist,char *brklist);
char   crack(int argc,char **argv,char *flags,int ign);


/* crack.c*/
extern int   arg_index;
extern char *arg_option;
extern char *pvcon;
