#include "real.h"
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

void bc4T(real ***t1,real ***t2, int i1,int i2,int j1,int j2,int k1,int k2,int bcw)
{

    //  Top and bottom faces for t //
    for (int row=j1; row<=j2; ++row) {
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
    for (int level=k1; level<=k2; ++level) {
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
    for (int level=k1; level<=k2; ++level) {
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

