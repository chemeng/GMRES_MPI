//
//  matmul.c
//  Askisi 1
//
//  Created by Tim on 4/4/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "methods_decl.h"

//Ypologismos Y=A*B matrix multiplication

void matmul(double **A, double **B, double **Y)
{
    double temp=0;
    int i=0,j=0,r=0;
    
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            temp=0;
            for (r=0; r<N; r++) {
                temp+=A[i][r]*B[r][j];
            }
            Y[i][j]=temp;
        }
    }
}