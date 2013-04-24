//
//  Krylov.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "parameters.h"

void Krylov(double **A,double *uj,double **Hm, double **u_base,int rank, 
        int nump,int local_N,int init,double *comm_tot)
{
    int  i=0,j=0,k=0;  
    double *w,*buffer,buff,help;
    double comm_time1=0,comm_time2=0;
    
    w=(double*)malloc(N*sizeof(double));
    buffer=(double*)malloc(local_N*sizeof(double));
    for (j=0;j<m;j++){
            //j einai to count g mas
            //vazei to kommati tou uj pou antistoixei sto process
        for(k=0;k<N;k++){
                uj[k]=u_base[k][j];
        }
        //matmul w=MATMUL(A,uj)
        for (i=0; i<local_N; i++) {
            w[i+init]=0;
            for (k=0; k<N; k++) {
                w[i+init]+=A[i][k]*uj[k];
            }
        }
        for (i=0;i<local_N; i++) {
            buffer[i]=w[i+init];
        }
        //broadcast tou w wste na to exoun oloi
        comm_time1=MPI_Wtime();
        MPI_Allgather(buffer, local_N, MPI_DOUBLE,w, local_N,
                      MPI_DOUBLE, MPI_COMM_WORLD);
        comm_time2=MPI_Wtime();
        *comm_tot += (comm_time2-comm_time1);
        
        for (i=0;i<=j;i++){
            for(k=0;k<N;k++){
                    uj[k]=u_base[k][i];
            }
            //DOT_PRODUCT(w,uj)
            buff=0;
            for (k=0; k<local_N; k++) {
                buff+=w[init+k]*uj[k+init];
            }
            comm_time1=MPI_Wtime();
            MPI_Allreduce(&buff,&help,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            comm_time2=MPI_Wtime();
            Hm[i][j]=help;
            *comm_tot+=(comm_time2-comm_time1);
            //end DOT product            
            for(k=0;k<=N;k++){
                    w[k]-=(Hm[i][j])*(uj[k]);
            }
        }
        buff=0;
        for (k=0;k<local_N;k++){
                buff+=pow(w[k+init],2);
        }
        comm_time1=MPI_Wtime();
        MPI_Allreduce(&buff,&help,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        comm_time2=MPI_Wtime();
        Hm[j+1][j]=help;
        *comm_tot+=(comm_time2-comm_time1);
        Hm[j+1][j]=sqrt(Hm[j+1][j]);
        if (Hm[j+1][j]==0) {
            break;
        }
        if(j<(m-1)){
            for (k=0;k<N;k++){
                    u_base[k][j+1]=w[k]/(Hm[j+1][j]);
            }
        } 
    }
}   

