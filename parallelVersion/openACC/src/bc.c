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

void bc(real ***restrict u,real ***restrict v,real ***restrict w, int i1,int i2,int j1,int j2,int k1,int k2,int bcw)
{
    
    #pragma acc parallel vector_length(VEC_LEN) present(u,v,w)
    {
    
        //  Top and bottom faces for u //
        #pragma acc loop gang collapse(3)
        for (int row=j1; row<=j2; ++row) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2+1; ++col) {
                //#pragma acc loop vector
                for (int k=1; k<=bcw; ++k) {        
                    u[k1-k][row][col] = u[k1][row][col];
                    u[k2+k][row][col] = u[k2][row][col];
                } // end for //
            } // end for //
        } // end for //

        //  East and West faces for u  //
        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2; ++row) {
                //#pragma acc loop vector
                for (int i=0; i<=bcw; ++i) {        
                    u[level][row][i1-i]   = -u[level][row][i1+1+i];
                    u[level][row][i2+1+i] = -u[level][row][i2-i];
                } // end for //
            } // end for //
        } // end for //

        //  North and South faces for u (periodic) //
        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2+1; ++col) {
                //#pragma acc loop vector
                for (int j=1; j<=bcw; ++j) {        
                    u[level][j1-j][col] = u[level][j2+1-j][col];
                    u[level][j2+j][col] = u[level][j1+j-1][col];
                } // end for //
            } // end for //
        } // end for //

    /////////////////////////////////////////////////////////////////////////////

        //  Top and bottom faces for v //
        #pragma acc loop gang collapse(3)
        for (int row=j1; row<=j2+1; ++row) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                //#pragma acc loop vector
                for (int k=1; k<=bcw; ++k) {        
                    v[k1-k][row][col] = v[k1][row][col];
                    v[k2+k][row][col] = v[k2][row][col];
                } // end for //
            } // end for //
        } // end for //


        //  East and West faces for v  //
        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2+1; ++row) {
                //#pragma acc loop vector
                for (int i=1; i<=bcw; ++i) {        
                    v[level][row][i1-i] = v[level][row][i1+i];
                    v[level][row][i2+i] = v[level][row][i2-i];
                } // end for //
            } // end for //
        } // end for //

        //  North and South faces for v (periodic) //
        #pragma acc loop gang collapse(2)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                v[level][j2+1][col] = v[level][j1][col];
                #pragma acc loop vector
                for (int j=1; j<=bcw; ++j) {        
                    v[level][j1-j][col]   = v[level][j2+1-j][col];
                    v[level][j2+1+j][col] = v[level][j1+j][col];
                } // end for //
            } // end for //
        } // end for //

    /////////////////////////////////////////////////////////////////////////////

        //  Top and bottom faces for w //
        #pragma acc loop gang collapse(3)
        for (int row=j1; row<=j2; ++row) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                //w[k1][row][col] = 0.0;
                //w[k2+1][row][col] = 0.0;
                //#pragma acc loop vector       
                for (int k=0; k<=bcw; ++k) {        
                    w[k1-k][row][col]   = 0.0;
                    w[k2+1+k][row][col] = 0.0;
                } // end for //
            } // end for //
        } // end for //



        //  East and West faces for w  //
        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2+1; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2; ++row) {
                //#pragma acc loop vector
                for (int i=1; i<=bcw; ++i) {        
                    w[level][row][i1-i] = w[level][row][i1+i];
                    w[level][row][i2+i] = w[level][row][i2-i];
                } // end for //
            } // end for //
        } // end for //


        //  North and South faces for w  (periodic) //
        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2+1; ++level) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                //#pragma acc loop vector
                for (int j=1; j<=bcw; ++j) {        
                    w[level][j1-j][col] = w[level][j2+1-j][col];
                    w[level][j2+j][col] = w[level][j1+j-1][col];
                } // end for //
            } // end for //
        } // end for //
    } // end of parallel region
} // end of bc() //

