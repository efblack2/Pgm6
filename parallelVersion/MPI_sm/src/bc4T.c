#include <mpi.h> 
#include "real.h"
#include "myMPI.h"

/*
 * ============================ bc =====================
 * BC sets the boundary conditions
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		main array size, not including
 *				extra 'ghost' zones/points
 */

void bc4T(real ***t1,real ***t2, int i1,int i2,int j1,int j2,int k1,int k2,int bcw, MPI_Comm sm_comm)
{
    int myRank, commSize;
    MPI_Comm_rank(sm_comm,&myRank);
    MPI_Comm_size(sm_comm,&commSize);


    const int jj1 = BLOCK_LOW (myRank,commSize,(1+j2-j1)) + bcw;
    const int jj2 = BLOCK_HIGH(myRank,commSize,(1+j2-j1)) + bcw;
    const int kk1 = BLOCK_LOW (myRank,commSize,(1+k2-k1)) + bcw;
    const int kk2 = BLOCK_HIGH(myRank,commSize,(1+k2-k1)) + bcw;
    //printf("myRank: %3d, --> jj1: %3d, jj2: %3d, j1: %3d, j2: %3d\n", myRank,jj1,jj2, j1,j2);
    //printf("myRank: %3d, --> kk1: %3d, kk2: %3d, k1: %3d, k2: %3d\n", myRank,kk1,kk2, k1,k2);

    //  Top and bottom faces for t //
    //#pragma omp for nowait
    for (int row=jj1; row<=jj2; ++row) {
        for (int col=i1; col<=i2; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                t1[k1-k][row][col] = t1[k1][row][col];
                t1[k2+k][row][col] = t1[k2][row][col];
                t2[k1-k][row][col] = t1[k1][row][col];
                t2[k2+k][row][col] = t1[k2][row][col];
            } // end for //
        } // end for //
    } // end for //


    //  East and West faces for t  //
    //#pragma omp for nowait
    for (int level=kk1; level<=kk2; ++level) {
        for (int row=j1; row<=j2; ++row) {
            for (int i=1; i<=bcw; ++i) {        
//                t1[level][row][i1-i] = t1[level][row][i1];
//                t1[level][row][i2+i] = t1[level][row][i2];
                t1[level][row][i1-i] = t1[level][row][i1+i];
                t1[level][row][i2+i] = t1[level][row][i2-i];
                t2[level][row][i1-i] = t1[level][row][i1+i];
                t2[level][row][i2+i] = t1[level][row][i2-i];
            } // end for //
        } // end for //
    } // end for //



    //  North and South faces for t  (periodic) //
    //#pragma omp for  // using the implicit barrirer to syncromize all thread before returning
    for (int level=kk1; level<=kk2; ++level) {
        for (int col=i1; col<=i2; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                t1[level][j1-j][col] = t1[level][j2+1-j][col];
                t1[level][j2+j][col] = t1[level][j1+j-1][col];
                t2[level][j1-j][col] = t1[level][j2+1-j][col];
                t2[level][j2+j][col] = t1[level][j1+j-1][col];
            } // end for //
        } // end for //
    } // end for //

} // end of bc4T() //
