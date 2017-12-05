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

void bc(real ***u,real ***v,real ***w, int i1,int i2,int j1,int j2,int k1,int k2,int bcw, int myRank, int commSize)
{

    const int jj1 = BLOCK_LOW (myRank,commSize,(1+j2-j1)) + bcw;
    const int jj2 = BLOCK_HIGH(myRank,commSize,(1+j2-j1)) + bcw;
    const int kk1 = BLOCK_LOW (myRank,commSize,(1+k2-k1)) + bcw;
    const int kk2 = BLOCK_HIGH(myRank,commSize,(1+k2-k1)) + bcw;
    //printf("myRank: %3d, --> jj1: %3d, jj2: %3d, j1: %3d, j2: %3d\n", myRank,jj1,jj2, j1,j2);
    //printf("myRank: %3d, --> kk1: %3d, kk2: %3d, k1: %3d, k2: %3d\n", myRank,kk1,kk2, k1,k2);
    
    //  Top and bottom faces for u //
    //#pragma omp for nowait
    for (int row=jj1; row<=jj2; ++row) {
        for (int col=i1; col<=i2+1; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                u[k1-k][row][col] = u[k1][row][col];
                u[k2+k][row][col] = u[k2][row][col];
            } // end for //
        } // end for //
    } // end for //

    //  East and West faces for u  //
    //#pragma omp for nowait
    for (int level=kk1; level<=kk2; ++level) {
        for (int row=j1; row<=j2; ++row) {
            for (int i=0; i<=bcw; ++i) {        
                u[level][row][i1-i]   = -u[level][row][i1+1+i];
                u[level][row][i2+1+i] = -u[level][row][i2-i];
            } // end for //
        } // end for //
    } // end for //

    //  North and South faces for u (periodic) //
    //#pragma omp for nowait
    for (int level=kk1; level<=kk2; ++level) {
        for (int col=i1; col<=i2+1; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                u[level][j1-j][col] = u[level][j2+1-j][col];
                u[level][j2+j][col] = u[level][j1+j-1][col];
            } // end for //
        } // end for //
    } // end for //

/////////////////////////////////////////////////////////////////////////////


    const int jj3 = BLOCK_LOW (myRank,commSize,(2+j2-j1)) + bcw;
    const int jj4 = BLOCK_HIGH(myRank,commSize,(2+j2-j1)) + bcw;
    //printf("myRank: %3d, --> jj3: %3d, jj4: %3d, j1: %3d, j2: %3d\n", myRank,jj3,jj4, j1,j2);

    //  Top and bottom faces for v //
    //#pragma omp for nowait
    for (int row=jj3; row<=jj4; ++row) {
        for (int col=i1; col<=i2; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                v[k1-k][row][col] = v[k1][row][col];
                v[k2+k][row][col] = v[k2][row][col];
            } // end for //
        } // end for //
    } // end for //


    //  East and West faces for v  //
    //#pragma omp for nowait
    for (int level=kk1; level<=kk2; ++level) {
        for (int row=j1; row<=j2+1; ++row) {
            for (int i=1; i<=bcw; ++i) {        
                v[level][row][i1-i] = v[level][row][i1+i];
                v[level][row][i2+i] = v[level][row][i2-i];
            } // end for //
        } // end for //
    } // end for //
 
    //  North and South faces for v (periodic) //
    //#pragma omp for nowait
    for (int level=kk1; level<=kk2; ++level) {
        for (int col=i1; col<=i2; ++col) {
            v[level][j2+1][col] = v[level][j1][col];
            for (int j=1; j<=bcw; ++j) {        
                v[level][j1-j][col]   = v[level][j2+1-j][col];
                v[level][j2+1+j][col] = v[level][j1+j][col];
            } // end for //
        } // end for //
    } // end for //

/////////////////////////////////////////////////////////////////////////////


    //  Top and bottom faces for w //
    //#pragma omp for nowait
    for (int row=jj1; row<=jj2; ++row) {
        for (int col=i1; col<=i2; ++col) {
            //w[k1][row][col] = 0.0;
            //w[k2+1][row][col] = 0.0;        
            for (int k=0; k<=bcw; ++k) {        
                w[k1-k][row][col]   = 0.0;
                w[k2+1+k][row][col] = 0.0;
            } // end for //
        } // end for //
    } // end for //


    const int kk3 = BLOCK_LOW (myRank,commSize,(2+k2-k1)) + bcw;
    const int kk4 = BLOCK_HIGH(myRank,commSize,(2+k2-k1)) + bcw;
    //printf("myRank: %3d, --> kk3: %3d, kk4: %3d, k1: %3d, k2: %3d\n", myRank,kk3,kk4, k1,k2);

    //  East and West faces for w  //
    //#pragma omp for nowait
    for (int level=kk3; level<=kk4; ++level) {
        for (int row=j1; row<=j2; ++row) {
            for (int i=1; i<=bcw; ++i) {        
                w[level][row][i1-i] = w[level][row][i1+i];
                w[level][row][i2+i] = w[level][row][i2-i];
            } // end for //
        } // end for //
    } // end for //


    //  North and South faces for w  (periodic) //
    //#pragma omp for  // using the implicit barrirer to syncromize all thread before returning
    for (int level=kk3; level<=kk4; ++level) {
        for (int col=i1; col<=i2; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                w[level][j1-j][col] = w[level][j2+1-j][col];
                w[level][j2+j][col] = w[level][j1+j-1][col];
            } // end for //
        } // end for //
    } // end for //

} // end of bc() //
