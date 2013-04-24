//
//  matvec.c
//  Askisi 1
//
//  Created by Tim on 4/4/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "methods_decl.h"

//Ypologismos y=A*x matrix vector multiplication

void matvec(double **A, double *x, double *y)
{
    double temp=0;
    int i=0,j=0;
    
    for (i=0; i<N; i++) {
        temp=0;
        for (j=0; j<N; j++) {
            temp+=A[i][j]*x[j];
        }
        y[i]=temp;
    }
}
