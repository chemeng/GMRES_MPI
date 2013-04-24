//
//  inv_power.c
//  Askisi 1
//
//  Created by Tim on 4/6/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include "methods_decl.h"
#include "parameters.h"

double inv_power(double **A)
{
    double eigenvalue=0,*q,*diff,*z,eps=0,temp=0,help=0,error=1;
    int i=0;
    q=(double*)(malloc(N*sizeof(double)));
    diff=(double*)(malloc(N*sizeof(double)));
    z=(double*)(malloc(N*sizeof(double)));
    eps=pow(10, -5);
    //initial guess
    for (i=0; i<N; i++) {
        q[i]=1;
    }
    eigenvalue=0;
    while (error > eps) {
        error=0;
        //find z_k=invA*q_k-1 via solution of system A*z_k=q_k-1
        GMRES_m(q,A,z);
        temp=norm(z);
        printf("norm=%lf\n",temp);
        for (i=0; i<N; i++) {
            help=q[i];
            q[i]=z[i]/temp;
            diff[i]=help-q[i];
            error+=pow(diff[i], 2);
        }
        error=sqrt(error);
    }
    for (i=0; i<N; i++) {
        eigenvalue+=z[i]*q[i];
    }
    return eigenvalue;
    
}