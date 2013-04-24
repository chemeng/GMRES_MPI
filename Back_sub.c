//
//  Back_sub.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include "parameters.h"

void Back_sub(double *y,double *g,double **Hm)
{
    int i=0,j=0,k=0;
    double sum=0,help[m+1];
  
    for (i=0;i<(m+1);i++){
        help[i]=g[i];
    }
    for (i=(m-1);i>=0;i--){
        sum=help[i];
        if (i<(m-1)){
            for (j=(i+1);j<m;j++){
                sum-=(Hm[i][j])*help[j];
            }
        }
        help[i]=sum/(Hm[i][i]);
    }
    for (i=0;i<m;i++){
        y[i]=help[i];
    }
}
