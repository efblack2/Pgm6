#include <mpi.h> 
#include "real.h"
#include "prototypes.h"
#include "myMPI.h"

void pgf(real ***restrict p3,real ***restrict u,real ***restrict v,real ***restrict w,real ***restrict t,real ***restrict p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int i1,int i2,int j1,int j2,int k1,int k2,real tstep,real g,real thetaBar, real cs2, int bcw, MPI_Comm sm_comm, MPI_Win *sm_win_u,MPI_Win *sm_win_v,MPI_Win *sm_win_w)  
{

    int myRank, commSize;
    MPI_Comm_rank(sm_comm,&myRank);
    MPI_Comm_size(sm_comm,&commSize);

    const int kk1 = BLOCK_LOW (myRank,commSize,(1+k2-k1)) + bcw;
    const int kk2 = BLOCK_HIGH(myRank,commSize,(1+k2-k1)) + bcw;
    //printf("myRank: %3d, --> kk1: %3d, kk2: %3d, k1: %3d, k2: %3d\n", myRank,kk1,kk2, k1,k2);

    //#pragma omp for nowait
    for(int level=kk1; level<=kk2; ++level) {
        real tempx = -tstep/(ro_u[level]*dx);
        for(int row=j1; row<=j2; ++row) {
            for(int col=i1+1; col<=i2; ++col) {
                u[level][row][col] += tempx * (p1[level-k1][row-j1][col-i1]-p1[level-k1][row-j1][col-i1-1]);
            } // end for //
        } // end for //
    } // end for //
    
    //#pragma omp for nowait
    for(int level=kk1; level<=kk2; ++level) {
        real tempy = -tstep/(ro_u[level]*dy);
        for(int row=j1+1; row<=j2; ++row) {
            for(int col=i1; col<=i2; ++col) {
                v[level][row][col] += tempy * (p1[level-k1][row-j1][col-i1]-p1[level-k1][row-j1-1][col-i1]);
            } // end for //
        } // end for //
    } // end for //

    const int kk3 = BLOCK_LOW (myRank,commSize,(k2-k1)) + bcw+1;
    const int kk4 = BLOCK_HIGH(myRank,commSize,(k2-k1)) + bcw+1;
    //printf("myRank: %3d, --> kk3: %3d, kk4: %3d, k1: %3d, k2: %3d\n", myRank,kk3,kk4, k1,k2);

    real temp1 = tstep*g/(2.0*thetaBar);
    //#pragma omp for // using the implicit barrier
    for(int level=kk3; level<=kk4; ++level) {
        real tempz = -tstep/(ro_w[level]*dz);
        for(int row=j1; row<=j2; ++row) {
            for(int col=i1; col<=i2; ++col) {
                w[level][row][col] += (tempz * (p1[level-k1][row-j1][col-i1]-p1[level-k1-1][row-j1][col-i1]) + temp1*(t[level][row][col] + t[level-1][row][col] - 2.0*thetaBar) ) ;
            } // end for //
        } // end for //
    } // end for //

    
    MPI_Win_sync(*sm_win_u);
    MPI_Win_sync(*sm_win_v);
    MPI_Win_sync(*sm_win_w);
    MPI_Barrier(sm_comm);

    bc(u,v,w,i1,i2,j1,j2,k1,k2,bcw,myRank,commSize);
    
    MPI_Win_sync(*sm_win_u);
    MPI_Win_sync(*sm_win_v);
    MPI_Win_sync(*sm_win_w);
    MPI_Barrier(sm_comm);

    temp1 = -tstep*cs2/dz;
    //#pragma omp for nowait
    for(int level=kk1; level<=kk2; ++level) {
        real tempx = -tstep*cs2*ro_u[level]/dx;
        real tempy = -tstep*cs2*ro_u[level]/dy;
        for(int row=j1; row<=j2; ++row) {
            for(int col=i1; col<=i2; ++col) {
                p3[level-k1][row-j1][col-i1] = p1[level-k1][row-j1][col-i1] + 
                                            ( tempx * (u[level][row][col+1] - u[level][row][col]) +  
                                              tempy * (v[level][row+1][col] - v[level][row][col]) +  
                                              temp1 * ( ro_w[level+1]*w[level+1][row][col] -  ro_w[level]*w[level][row][col] ));
            } // end for //
        } // end for //
    } // end for //    
//    MPI_Barrier(sm_comm);
} // end of pgf() //
