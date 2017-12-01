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
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */

void advec1(real *q2,const real *q1,real *uv,real dx,real dt,int iW,int iE, int nxydim, char advection_type);

void advection(real ***__restrict__  q2,real ***u,real ***v,real ***w, real ***q1, real ***__restrict__  uu,real ***__restrict__  vv, real ***__restrict__  ww,real dx,real dy,real dz,real dt,int iW,int iE,int jS,int jN,int kB,int kT,int nxdim,int nydim,int nzdim,real tstep, char advection_type)


{
    real *__restrict__  fld;

    // advection in x //
    fld = (real *) malloc(nxdim*sizeof(real));
    
    for (int level=kB; level<=kT; ++level) {	
        for (int row=jS; row<=jN; ++row) {	
            if (advection_type != 'p') {
                for (int col=iW; col<=iE; ++col) {
                    fld[col] = 0.5*(u[level][row][col] + u[level][row][col+1]);  
                } // end for //
            } else {
                for (int col=iW; col<=iE+1; ++col) {
                    fld[col] = u[level][row][col];
                } // end for //
            } // end if //
            advec1(q2[level][row],q1[level][row],fld,dx,dt,iW,iE,nxdim,advection_type);
        } // end for //
    } // end for //
    free(fld);
    
    // advection in y //
    
    real *__restrict__  qIn,*__restrict__  qOut;

    qIn  = (real *) malloc(nydim*sizeof(real));
    qOut = (real *) malloc(nydim*sizeof(real));
    fld  = (real *) malloc(nydim*sizeof(real));
    
    for (int level=kB; level<=kT; ++level) {	
        for (int col=iW; col<=iE; ++col) {
            for (int row=0; row<nydim; ++row) {	
                qIn[row] = q2[level][row][col];
            } // end for //
            if (advection_type != 'p') {
                for (int row=jS; row<=jN; ++row) {
                    fld[row] = 0.5*(v[level][row][col] + v[level][row+1][col]);  
                } // end for //
            } else {
                for (int row=jS; row<=jN+1; ++row) {
                    fld[row] = v[level][row][col];
                } // end for //
            } // end if //
            advec1(qOut,qIn,fld,dy,dt,jS,jN,nydim,advection_type);
            for (int row=jS; row<=jN; ++row) {
                q2[level][row][col] = qOut[row];
            } // end for //
        } // end for //
    } // end for //
    free(qIn);
    free(qOut);
    free(fld);

    
     // advection in z //
    qIn  = (real *) malloc(nzdim*sizeof(real));
    qOut = (real *) malloc(nzdim*sizeof(real));
    fld  = (real *) malloc(nzdim*sizeof(real));
    
    for (int row=jS; row<=jN; ++row) {	
        for (int col=iW; col<=iE; ++col) {
            for (int level=0; level<nzdim; ++level) {	
                qIn[level] = q2[level][row][col];
            } // end for //
            if (advection_type != 'p') {
                for (int level=kB; level<=kT; ++level) {
                    fld[level] = 0.5*(w[level][row][col] + w[level+1][row][col]);  
                } // end for //
            } else {
                for (int level=kB; level<=kT+1; ++level) {
                    fld[level] = w[level][row][col];
                } // end for //
            } // end if //
            advec1(qOut,qIn,fld,dz,dt,kB,kT,nzdim,advection_type);
            for (int level=kB; level<=kT; ++level) {
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
    
    for (int level=kB; level<=kT; ++level) {	
        for (int row=jS; row<=jN; ++row) {	
            for (int col=iW; col<=iE+1; ++col) {	
                u[level][row][col] += tempx *((uu[level][row][col]+uu[level][row][col-1])     * (uu[level][row][col]-uu[level][row][col-1])   +
                                              (uu[level][row][col+1]+uu[level][row][col])     * (uu[level][row][col+1]-uu[level][row][col]))  +
                                      tempy *((vv[level][row][col]+vv[level][row][col-1])     * (uu[level][row][col]-uu[level][row-1][col])   +   
                                              (vv[level][row+1][col]+vv[level][row+1][col-1]) * (uu[level][row+1][col]-uu[level][row][col]))  +                
                                      tempz *((ww[level][row][col]+ww[level][row][col-1])     * (uu[level][row][col]-uu[level-1][row][col])   +   
                                              (ww[level+1][row][col]+ww[level+1][row][col-1]) * (uu[level+1][row][col]-uu[level][row][col]));
            } // end for //
        } // end for //
    } // end for //
    
    for (int level=kB; level<=kT; ++level) {	
        for (int row=jS; row<=jN+1; ++row) {	
            for (int col=iW; col<=iE; ++col) {	            
                v[level][row][col] += tempx *((uu[level][row][col]+uu[level][row-1][col])     * (vv[level][row][col]-vv[level][row][col-1])   + 
                                              (uu[level][row][col+1]+uu[level][row-1][col+1]) * (vv[level][row][col+1]-vv[level][row][col]))  + 
                                      tempy *((vv[level][row][col]+vv[level][row-1][col])     * (vv[level][row][col]-vv[level][row-1][col])   + 
                                              (vv[level][row+1][col]+vv[level][row][col])     * (vv[level][row+1][col]-vv[level][row][col]))  +               
                                      tempz *((ww[level][row][col]+ww[level][row-1][col])     * (vv[level][row][col]-vv[level-1][row][col])   + 
                                              (ww[level+1][row][col]+ww[level+1][row-1][col]) * (vv[level+1][row][col]-vv[level][row][col]));
            } // end for //
        } // end for //
    } // end for //

    for (int level=kB; level<=kT+1; ++level) {	
        for (int row=jS; row<=jN; ++row) {	
            for (int col=iW; col<=iE; ++col) {
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
