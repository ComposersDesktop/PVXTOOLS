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
/*
 * sfexpr - similar to expr() routine from libfrm.a, but optimized for
 * parsing command line for csound.  No major differences, some
 * features left out, such as sin(), ln(), etc.  Postoperators added:
 * s - seconds, S - samples, K - 1024, ms - milliseconds, m - minutes,
 * dB. 
 */
/* modified for doubles */

/*------------------------------------------------------------------------------

        Sfexpr: Updated by T.Wishart and N.Laviers 1991.

------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cmdutils.h"

int
sffield(char **input,char *string,char *iglist,char *brklist);

double 
sfexpr(char *string,double timefac)  
{
        char *rpnp, rpn[LRPN], *it, item[40], sym[40];
        int nops, n, top = 0;
        double stack[LRPN];

        rpnp = rpn;

        strncpy(rpn,polish(string,UNOPS,BINOPS,SFPOSTOPS),LRPN);

        while(strlen(rpnp)){
        
        sffield(&rpnp, item, "", ",");
        it=item;
        sffield(&it, sym, "", "$");
        nops = atoi(it);
   
        if (!nops){   /* it is a constant number */

                if(!strchr(sym,'.') && *sym == '0'){ 
                        /* no '.' and '0' prefix */

                    if( *(sym+1) == 'x' || *(sym+1) == 'X' ) 
                    sscanf(sym+2, "%x", &n);
                else 
                    sscanf(sym, "%o", &n);

                    stack[++top] = n;
                } 
            else 
                sscanf(sym,"%lf",&stack[++top]);
                stack[top] *= timefac;
                continue;
            }

        /* from here on, deal with operators only */

        if (!strcmp(sym,"-") && nops == 1){
                stack[top] = -stack[top];
            continue;
        }else if (nops == 2){

                switch(sym[0]){
                    case '-':
                    stack[top-1] = stack[top-1] - stack[top];
                    break;

                    case '+':
                    stack[top-1] = stack[top-1] + stack[top];
                    break;

                    case '*':
                    stack[top-1] = stack[top-1] * stack[top];
                    break;
                
                case '/':
                    stack[top-1] = stack[top-1] / stack[top];
                    break;

                    case '^':
                    stack[top-1] = pow(stack[top-1], stack[top]);
                    break;
                
                case '%':
                    stack[top-1] = (int) stack[top-1] % (int) stack[top];
            }
                top--;
        }else if (nops == 1){
            
            if(!strcmp(sym,"dB")){
                    stack[top] /= timefac;
                    stack[top] = pow( (double) 10.0, stack[top]/20.);
                    continue;
            }

                if(!strcmp(sym,"K")){
                    stack[top] *= 1024.0;
                    continue;
                }

                if(!strcmp(sym,"k")){
                    stack[top] *= 1000.0;
                continue;
            }

                if(!strcmp(sym,"S")){
                    stack[top] /= timefac;
                    continue;
            }

                if(!strcmp(sym,"ms")){
                    stack[top] *= 0.001;
                continue;
            }

                if(!strcmp(sym,"m")){
                    stack[top] *= 60.0;
                continue;
                }

                if(!strcmp(sym,"s")){
                    continue;
                }

            }
    }
        return(stack[top]);
}



int
sffield(char **input,char *string,char *iglist,char *brklist) 
{

    int leading = 1, c;

    while ((c = *(*input)++) != 0) {

            if (leading && c==' ')
            continue;

            if (strchr(iglist,c))
                continue;

            if (!strchr(brklist,c)){
                *string++ = c; 
                leading = 0; 
                continue;
            }else{
            *string = 0; 
                return(0);
        }
    }
        return(0); 
}

