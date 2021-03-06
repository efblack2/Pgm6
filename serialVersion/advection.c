#include <stdio.h>
#include <stdlib.h>
#include "real.h"

/*
 * ======================= advection ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
 *	c	real		true speed of wave
 *	dx	real		grid spacing
 *	dt	real		time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */

void advec1(real *q2,const real *q1,real *uv,real dx,real dt,int i1,int i2, char advection_type);

void advection(real ***restrict q2,real ***u,real ***v,real ***w, real ***q1, real ***restrict uu,real ***restrict vv, real ***restrict ww,real dx,real dy,real dz,real dt,int i1,int i2,int j1,int j2,int k1,int k2,int nxdim,int nydim,int nzdim,real tstep, char advection_type)


{
    real *restrict fld;

    // advection in x //
    fld =  malloc(nxdim*sizeof(real));
    
    for (int level=k1; level<=k2; ++level) {	
        for (int row=j1; row<=j2; ++row) {	
            if (advection_type != 'p') {
                for (int col=i1; col<=i2; ++col) {
                    fld[col] = 0.5*(u[level][row][col] + u[level][row][col+1]);  
                } // end for //
            } else {
                for (int col=i1; col<=i2+1; ++col) {
                    fld[col] = u[level][row][col];
                } // end for //
            } // end if //
            advec1(q2[level][row],q1[level][row],fld,dx,dt,i1,i2,advection_type);
        } // end for //
    } // end for //
    free(fld);
    
    // advection in y //
    
    real *restrict qIn,*restrict qOut;

    qIn =  malloc(nydim*sizeof(real));
    qOut = malloc(nydim*sizeof(real));
    fld =  malloc(nydim*sizeof(real));
    
    for (int level=k1; level<=k2; ++level) {	
        for (int col=i1; col<=i2; ++col) {
            for (int row=0; row<nydim; ++row) {	
                qIn[row] = q2[level][row][col];
            } // end for //
            if (advection_type != 'p') {
                for (int row=j1; row<=j2; ++row) {
                    fld[row] = 0.5*(v[level][row][col] + v[level][row+1][col]);  
                } // end for //
            } else {
                for (int row=j1; row<=j2+1; ++row) {
                    fld[row] = v[level][row][col];
                } // end for //
            } // end if //
            advec1(qOut,qIn,fld,dy,dt,j1,j2,advection_type);
            for (int row=j1; row<=j2; ++row) {
                q2[level][row][col] = qOut[row];
            } // end for //
        } // end for //
    } // end for //
    free(qIn);
    free(qOut);
    free(fld);

    
     // advection in z //
    qIn =  malloc(nzdim*sizeof(real));
    qOut = malloc(nzdim*sizeof(real));
    fld =  malloc(nzdim*sizeof(real));
    
    for (int row=j1; row<=j2; ++row) {	
        for (int col=i1; col<=i2; ++col) {
            for (int level=0; level<nzdim; ++level) {	
                qIn[level] = q2[level][row][col];
            } // end for //
            if (advection_type != 'p') {
                for (int level=k1; level<=k2; ++level) {
                    fld[level] = 0.5*(w[level][row][col] + w[level+1][row][col]);  
                } // end for //
            } else {
                for (int level=k1; level<=k2+1; ++level) {
                    fld[level] = w[level][row][col];
                } // end for //
            } // end if //
            advec1(qOut,qIn,fld,dz,dt,k1,k2,advection_type);
            for (int level=k1; level<=k2; ++level) {
                q2[level][row][col] = qOut[level];
            } // end for //
        } // end for //
    } // end for //

    
    free(qIn);
    free(qOut);
    free(fld);

   	real tempx,tempy,tempz;
    tempx = -tstep/(4*dx);
    tempy = -tstep/(4*dy);
    tempz = -tstep/(4*dz);
    
    for (int level=k1; level<=k2; ++level) {	
        for (int row=j1; row<=j2; ++row) {	
            for (int col=i1; col<=i2+1; ++col) {	
                u[level][row][col] += tempx *((uu[level][row][col]+uu[level][row][col-1])     * (uu[level][row][col]-uu[level][row][col-1])   +
                                              (uu[level][row][col+1]+uu[level][row][col])     * (uu[level][row][col+1]-uu[level][row][col]))  +
                                      tempy *((vv[level][row][col]+vv[level][row][col-1])     * (uu[level][row][col]-uu[level][row-1][col])   +   
                                              (vv[level][row+1][col]+vv[level][row+1][col-1]) * (uu[level][row+1][col]-uu[level][row][col]))  +                
                                      tempz *((ww[level][row][col]+ww[level][row][col-1])     * (uu[level][row][col]-uu[level-1][row][col])   +   
                                              (ww[level+1][row][col]+ww[level+1][row][col-1]) * (uu[level+1][row][col]-uu[level][row][col]));
            } // end for //
        } // end for //
    } // end for //
    
    for (int level=k1; level<=k2; ++level) {	
        for (int row=j1; row<=j2+1; ++row) {	
            for (int col=i1; col<=i2; ++col) {	            
                v[level][row][col] += tempx *((uu[level][row][col]+uu[level][row-1][col])     * (vv[level][row][col]-vv[level][row][col-1])   + 
                                              (uu[level][row][col+1]+uu[level][row-1][col+1]) * (vv[level][row][col+1]-vv[level][row][col]))  + 
                                      tempy *((vv[level][row][col]+vv[level][row-1][col])     * (vv[level][row][col]-vv[level][row-1][col])   + 
                                              (vv[level][row+1][col]+vv[level][row][col])     * (vv[level][row+1][col]-vv[level][row][col]))  +               
                                      tempz *((ww[level][row][col]+ww[level][row-1][col])     * (vv[level][row][col]-vv[level-1][row][col])   + 
                                              (ww[level+1][row][col]+ww[level+1][row-1][col]) * (vv[level+1][row][col]-vv[level][row][col]));
            } // end for //
        } // end for //
    } // end for //

    for (int level=k1; level<=k2+1; ++level) {	
        for (int row=j1; row<=j2; ++row) {	
            for (int col=i1; col<=i2; ++col) {
                w[level][row][col] += tempx *((uu[level][row][col]+uu[level-1][row][col])     * (ww[level][row][col]-ww[level][row][col-1])  + 
                                              (uu[level][row][col+1]+uu[level-1][row][col+1]) * (ww[level][row][col+1]-ww[level][row][col])) +
                                      tempy *((vv[level][row][col]+vv[level-1][row][col])     * (ww[level][row][col]-ww[level][row-1][col])  + 
                                              (vv[level][row+1][col]+vv[level-1][row+1][col]) * (ww[level][row+1][col]-ww[level][row][col])) +
                                      tempz *((ww[level][row][col]+ww[level-1][row][col])     * (ww[level][row][col]-ww[level-1][row][col])  + 
                                              (ww[level+1][row][col]+ww[level][row][col])     * (ww[level+1][row][col]-ww[level][row][col]));
            } // end for //
        } // end for //
    } // end for //
} // end of advection() //
