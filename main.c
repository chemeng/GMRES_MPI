//
//  main.c
//  Homework 2
//
//  Created by Tim on 4/2/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h>
#include "methods_decl.h"
#include "parameters.h"
#include "inv_power.c"
#include "norm.c"
#include "matmul.c"
#include "matvec.c"
#include "matrix_fill.c"
#include "GMRES_m.c"
#include "Krylov.c"
#include "Back_sub.c"
#include "Least_sq.c"


int main (int argc,char* argv[])
{
    double *x,*y,**A,**S,**invS,**D,**T,temp=0,help=0,eigen=0,**local_A;
    int i=0,j=0,local_N,temp_div,temp_mod,init,fin;
    int nump=0,rank=0;
    clock_t start1, end1 ;
    double time1,time2;
    int mem1=0,mem2=0;
    //MPI world build
    MPI_Init( &argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nump);
    //allocation
    x=(double*)(malloc(N*sizeof(double)));
    y=(double*)(malloc(N*sizeof(double)));
    //allocation of 2D A[N][N]
    A=(double**)(malloc(N*sizeof(double*)));
    S=(double**)(malloc(N*sizeof(double*)));
    D=(double**)(malloc(N*sizeof(double*)));
    invS=(double**)(malloc(N*sizeof(double*)));
    T=(double**)(malloc(N*sizeof(double*)));
    for (i=0; i<N; i++) {
        S[i]=(double*)(malloc(N*sizeof(double)));
        D[i]=(double*)(malloc(N*sizeof(double)));
        invS[i]=(double*)(malloc(N*sizeof(double)));
        T[i]=(double*)(malloc(N*sizeof(double)));
    } 
    /*
     #####  Assemble of matrix A   #####
    */
    //fill matrices S,D,invS
    matrix_fill(x,S,D,invS);
    
    //matrix matrix multiplication T=S*D
    matmul(S,D,T);
    free(D);
    free(S);
    //matrix matrix multiplication A=T*invS
    for (i=0; i<N; i++) {
        A[i]=(double*)(malloc(N*sizeof(double)));
    }
    matmul(T,invS,A);
    //matrix vector multiplication y=Ax
    matvec(A,x,y);
    free(T);
    free(invS);
    
    //Partition of matrix A by lines
    //All processes have vector b
    
    //Distribution of load (dimension N_local) 
    temp_mod=N%nump;
    temp_div=N/nump;
    if (temp_mod==0){
        init=rank*temp_div;
        fin=(rank+1)*temp_div-1;
    }
    else if (temp_mod!=0){
        if (rank!=(nump-1)){
            init=(rank*temp_div+rank);
            fin=((rank+1)*(temp_div+1)-1);
        }
        else if (rank==(nump-1)){
            init=(rank*temp_div+rank);
            fin=(N-1);
        }
    }
    local_N=fin-init+1;
    //#####  GMRES   #####
    double *b;
    b=(double*)(malloc(N*sizeof(double)));
    for (i=0; i<N; i++) {
        x[i]=0;
        b[i]=1;
    }
    //allocation of 2D local array A
    local_A=(double**)(malloc(local_N*sizeof(double*)));
    for (i=0; i<local_N; i++) {
        local_A[i]=(double*)(malloc(N*sizeof(double)));
    }
    //Fill local arrays with corresponding data chunks
    for (i=0; i<local_N; i++) {
        for (j=0; j<N; j++) {
            local_A[i][j]=A[i+init][j];
        }
    }
    free(A);
    //Solve Ax=b with parallel GMRES
    GMRES_m(b,local_A,x,rank,nump,local_N,init,fin);
    if (rank==0) {
        printf("GMRES(m) DONE, m=%d\n",m);
        printf("Result printed by Rank %d: x[N/2]=%.10lf\n",rank,x[N/2]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}


















