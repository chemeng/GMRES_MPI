//
//  matrix_fill.c
//  Askisi 1
//
//  Created by Tim on 4/4/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parameters.h"
#include "methods_decl.h"

void matrix_fill(double *x, double **S, double **D, double **invS)
{
    int i=0,j=0;
    double temp;
    
    for (i=0; i<N; i++) {
        x[i]=i+1;
        for (j=0; j<N; j++) {
            S[i][j]=0;
            D[i][j]=0;
            invS[i][j]=0;
            if (i==j) {
                S[i][j]=1;
                D[i][j]=i+1;
                invS[i][j]=1;
            }
            if ((j-i)==1) {
                S[i][j]=beta;
            }
            if (i<=j) {
                temp=(double) j-i;
                invS[i][j]=pow(-beta, temp);
            }
        }
    }
    
}