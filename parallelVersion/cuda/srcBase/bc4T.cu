#include "real.h"
/*
 * ============================ bc =====================
 * BC sets the boundary conditions
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current time step
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		main array size, not including
 *				extra 'ghost' zones/points
 */

void bc4T(real ***t1,real ***t2, int iW,int iE,int jS,int jN,int kB,int kT,int bcw)
{

    //  Top and bottom faces for t //
    for (int row=jS; row<=jN; ++row) {
        for (int col=iW; col<=iE; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                t1[kB-k][row][col] = t1[kB][row][col];
                t1[kT+k][row][col] = t1[kT][row][col];
                t2[kB-k][row][col] = t1[kB][row][col];
                t2[kT+k][row][col] = t1[kT][row][col];                
            } // end for //
        } // end for //
    } // end for //


    //  East and West faces for t  //
    for (int level=kB; level<=kT; ++level) {
        for (int row=jS; row<=jN; ++row) {
            for (int i=1; i<=bcw; ++i) {        
//                t1[level][row][iW-i] = t1[level][row][iW];
//                t1[level][row][iE+i] = t1[level][row][iE];
                t1[level][row][iW-i] = t1[level][row][iW+i];
                t1[level][row][iE+i] = t1[level][row][iE-i];
                t2[level][row][iW-i] = t1[level][row][iW+i];
                t2[level][row][iE+i] = t1[level][row][iE-i];
            } // end for //
        } // end for //
    } // end for //



    //  North and South faces for t  (periodic) //
    for (int level=kB; level<=kT; ++level) {
        for (int col=iW; col<=iE; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                t1[level][jS-j][col] = t1[level][jN+1-j][col];
                t1[level][jN+j][col] = t1[level][jS+j-1][col];
                t2[level][jS-j][col] = t1[level][jN+1-j][col];
                t2[level][jN+j][col] = t1[level][jS+j-1][col];
            } // end for //
        } // end for //
    } // end for //
} // end of bc4T() //

