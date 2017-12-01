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

void bc(real ***u,real ***v,real ***w, int iW,int iE,int jS,int jN,int kB,int kT,int bcw)
{
    //  Top and bottom faces for u //
    for (int row=jS; row<=jN; ++row) {
        for (int col=iW; col<=iE+1; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                u[kB-k][row][col] = u[kB][row][col];
                u[kT+k][row][col] = u[kT][row][col];
            } // end for //
        } // end for //
    } // end for //

    //  East and West faces for u  //
    for (int level=kB; level<=kT; ++level) {
        for (int row=jS; row<=jN; ++row) {
            for (int i=0; i<=bcw; ++i) {        
                u[level][row][iW-i]   = -u[level][row][iW+1+i];
                u[level][row][iE+1+i] = -u[level][row][iE-i];
            } // end for //
        } // end for //
    } // end for //

    //  North and South faces for u (periodic) //
    for (int level=kB; level<=kT; ++level) {
        for (int col=iW; col<=iE+1; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                u[level][jS-j][col] = u[level][jN+1-j][col];
                u[level][jN+j][col] = u[level][jS+j-1][col];
            } // end for //
        } // end for //
    } // end for //

/////////////////////////////////////////////////////////////////////////////


    //  Top and bottom faces for v //
    for (int row=jS; row<=jN+1; ++row) {
        for (int col=iW; col<=iE; ++col) {
            for (int k=1; k<=bcw; ++k) {        
                v[kB-k][row][col] = v[kB][row][col];
                v[kT+k][row][col] = v[kT][row][col];
            } // end for //
        } // end for //
    } // end for //


    //  East and West faces for v  //
    for (int level=kB; level<=kT; ++level) {
        for (int row=jS; row<=jN+1; ++row) {
            for (int i=1; i<=bcw; ++i) {        
                v[level][row][iW-i] = v[level][row][iW+i];
                v[level][row][iE+i] = v[level][row][iE-i];
            } // end for //
        } // end for //
    } // end for //

    //  North and South faces for v (periodic) //
    for (int level=kB; level<=kT; ++level) {
        for (int col=iW; col<=iE; ++col) {
            v[level][jN+1][col] = v[level][jS][col];
            for (int j=1; j<=bcw; ++j) {        
                v[level][jS-j][col]   = v[level][jN+1-j][col];
                v[level][jN+1+j][col] = v[level][jS+j][col];
            } // end for //
        } // end for //
    } // end for //

/////////////////////////////////////////////////////////////////////////////


    //  Top and bottom faces for w //
    for (int row=jS; row<=jN; ++row) {
        for (int col=iW; col<=iE; ++col) {
            //w[kB][row][col] = 0.0;
            //w[kT+1][row][col] = 0.0;        
            for (int k=0; k<=bcw; ++k) {        
                w[kB-k][row][col]   = 0.0;
                w[kT+1+k][row][col] = 0.0;
            } // end for //
        } // end for //
    } // end for //



    //  East and West faces for w  //
    for (int level=kB; level<=kT+1; ++level) {
        for (int row=jS; row<=jN; ++row) {
            for (int i=1; i<=bcw; ++i) {        
                w[level][row][iW-i] = w[level][row][iW+i];
                w[level][row][iE+i] = w[level][row][iE-i];
            } // end for //
        } // end for //
    } // end for //


    //  North and South faces for w  (periodic) //
    for (int level=kB; level<=kT+1; ++level) {
        for (int col=iW; col<=iE; ++col) {
            for (int j=1; j<=bcw; ++j) {        
                w[level][jS-j][col] = w[level][jN+1-j][col];
                w[level][jN+j][col] = w[level][jS+j-1][col];
            } // end for //
        } // end for //
    } // end for //

} // end of bc() //
