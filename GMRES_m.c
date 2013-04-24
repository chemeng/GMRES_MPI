//
//  gmres.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "parameters.h"
#include "methods_decl.h"

//int m=0;

void GMRES_m(double *b,double **A,double *x,int rank, int nump,int local_N,int init,int fin)
{
    int i=0,j=0,k=0,iter=0;
    double vita=0,eps=0.000000001,sum=0,buff,*buffer;
    double e[m+1],y[m],g[m+1];
    double *x0,*r0,*uj,*temp,norm=1,norm0=1;
    double **Hm,**u_base;
    double tot_time=0,init_time=0,comm_time1=0,comm_time2=0,comm_tot=0;
    //allocation wste na apothikeutoun sto heap
    init_time=MPI_Wtime();
    x0=(double*)calloc(N,sizeof(double));
    r0=(double*)malloc(N*sizeof(double));
    buffer=(double*)malloc(local_N*sizeof(double));
    uj=(double*)malloc(N*sizeof(double));
    temp=(double*)malloc(N*sizeof(double));    
    //allocation of 2D A[N][M]
    Hm=(double**)(malloc((m+1)*sizeof(double*)));
    for (i=0; i<(m+1); i++) {
        Hm[i]=(double*)(malloc(m*sizeof(double)));
    }
    u_base=(double**)(malloc(N*sizeof(double*)));
    for (i=0; i<N; i++) {
        x[i]=0;
        u_base[i]=(double*)(malloc(m*sizeof(double)));
    }
    iter=1;
    while ((norm/norm0)>=eps) {
        for (i=0; i<N; i++) {
            x0[i]=x[i];
            uj[i]=0;
            r0[i]=0;
            temp[i]=0;
        }
        //Initialization
        for (i=0; i<N; i++) {
            for (j=0; j<m; j++) {
                u_base[i][j]=0;
            }
        }
        for (i=0; i<(m+1); i++) {
            e[i]=0;
            g[i]=0;
            for (j=0; j<m; j++) {
                Hm[i][j]=0;
                y[j]=0;
            }
        }
        e[0]=1;
        //r0=b-A*x0 
        for (i=0; i<local_N; i++) {
            sum=0;
            for (j=0; j<N; j++) {
                sum+=A[i][j]*x0[j];
            }
            r0[i+init]=b[i+init]-sum;
        }                 
        vita=0;
        buff=0;
            //norm r0
        for (i=0; i<local_N; i++) {
                vita+=pow(r0[init+i],2);  
        }
        comm_time1=MPI_Wtime();
        MPI_Allreduce(&vita,&buff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        comm_time2=MPI_Wtime();
        comm_tot+=(comm_time2-comm_time1);
        
        vita=sqrt(buff);
        if (iter==1) {
            norm0=vita;
        }
        norm=vita;
            //Initial vector u_j
        for (i=0;i<local_N; i++) {
            uj[i+init]=r0[i+init]/vita;
            buffer[i]=uj[i+init];
        }
        //broadcast uj to all processes
        comm_time1=MPI_Wtime();
        MPI_Allgather(buffer, local_N, MPI_DOUBLE,uj, local_N,
                    MPI_DOUBLE, MPI_COMM_WORLD);
        comm_time2=MPI_Wtime();
        comm_tot+=(comm_time2-comm_time1);
        //Save u1 to Krylov matrix
        for (i=0; i<N; i++) {
            u_base[i][0]=uj[i];
        }
            //vector g=b*e1
        for (i=0; i<(m+1); i++) {
            g[i]=vita*e[i];
        }
        //Build Krylov subspace and solve least squares
        Krylov(A,uj,Hm,u_base,rank,nump,local_N,init,&comm_tot);
        Least_sq(g,Hm);
        Back_sub(y,g,Hm);
        //x=x0+MATMUL(u_base(N,m),y(m,1)))
        for (i=0; i<local_N; i++) {
            for (k=0; k<m; k++) {
                temp[i]+=u_base[i+init][k]*y[k];
            }
        }
        for (i=0; i<local_N; i++) {
            x[i+init]=x0[i+init]+temp[i];
        }
        for (i=0;i<local_N; i++) {
            buffer[i]=x[i+init];
        }
        //Brodcast u_j to all processes
        comm_time1=MPI_Wtime();
        MPI_Allgather(buffer, local_N, MPI_DOUBLE,x, local_N,
                      MPI_DOUBLE, MPI_COMM_WORLD);
        comm_time2=MPI_Wtime();
        comm_tot+=(comm_time2-comm_time1);
        comm_time2=comm_tot;

        iter++;
    }
    free(x0);
    free(r0);
    free(uj);
    free(temp);
    free(Hm);
    free(u_base);
    comm_time1=MPI_Wtime()-init_time;
    MPI_Allreduce(&comm_time1,&tot_time,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&comm_time2,&comm_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (rank==0) {
        printf("N=%d  NUMP=%d\n",N,nump);
        printf("Total time=%.10lf\n",tot_time/nump);
        printf("Communication time=%.10lf\n",comm_tot/nump);
        printf("Computation time=%.10lf\n",(tot_time-comm_tot)/nump);
    }
   // }
}
                                                          
                                                          