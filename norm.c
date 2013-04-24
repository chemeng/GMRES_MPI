//
//  norm.c
//  Askisi 1
//
//  Created by Tim on 4/6/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "parameters.h"
#include "methods_decl.h"

double norm(double *x)
{
    int i=0;
    double temp,res=0;
    temp=0;
    for (i=0; i<N; i++) {
        temp+=pow(x[i], 2);
    }
    res=sqrt(temp);
    return res;
}