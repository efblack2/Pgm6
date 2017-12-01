#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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


#pragma acc routine vector
void advec1(real *q2,const real *q1,real *uv,real *flux, real dx,real dt,int i1,int i2, int nxydim, char advection_type);

void advection(real ***restrict q2, real ***restrict u ,real ***restrict v ,real ***restrict w ,
               real ***restrict q1, real ***restrict uu,real ***restrict vv,real ***restrict ww,
               real dx,real dy,real dz,real dt,int i1,int i2,int j1,int j2,int k1,int k2,int nxdim,int nydim,int nzdim,real tstep, char advection_type)
{
    real *restrict fld;
    real *restrict flux;

    // advection in x //
    fld  =  malloc(nxdim*sizeof(real));
    flux =  malloc(nxdim*sizeof(real));

    #pragma acc parallel vector_length(VEC_LEN) present(q2,u,q1) private (fld[0:nxdim], flux[0:nxdim])
    {
        #pragma acc loop gang collapse(2)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2; ++row) {
                if (advection_type != 'p') {
                    #pragma acc loop vector
                    for (int col=i1; col<=i2; ++col) {
                        fld[col] = 0.5*(u[level][row][col] + u[level][row][col+1]);
                    } // end for //
                } else {
                    #pragma acc loop vector
                    for (int col=i1; col<=i2+1; ++col) {
                        fld[col] = u[level][row][col];

                        real r = fabs(fld[col]*dt/dx);
                        if (fld[col] >=0.0) {
                            flux[col]= r*(q1[level][row][col-1] + 0.25*(1-r)*(q1[level][row][col]-q1[level][row][col-2]));
                        } else {
                            flux[col]= r*(-q1[level][row][col]  + 0.25*(1-r)*(q1[level][row][col+1] - q1[level][row][col-1]));
                        } // end if //
                    } // end for //
                    #pragma acc loop vector
                    for (int col=i1; col<=i2; ++col) {
                        q2[level][row][col] = q1[level][row][col] * (1.0 + (fld[col+1] - fld[col])*dt/dx ) - (flux[col+1] - flux[col]);
                        //q2[level][row][col] = q1[level][row][col] - (flux[col+1] - flux[col]) + dt/dx*q1[level][row][col]*(fld[col+1] - fld[col]);
                    } // end for //
                } // end if //
                //advec1(q2[level][row],q1[level][row],fld,flux,dx,dt,i1,i2,nxdim,advection_type);
            } // end for //
        } // end for //
    } // end of parallel region //
    free(fld);
    free(flux);


    real *restrict qIn;

    qIn  = malloc(nydim*sizeof(real));
    fld  = malloc(nydim*sizeof(real));
    flux = malloc(nydim*sizeof(real));
    #pragma acc parallel vector_length(VEC_LEN) present(q2,v,q1) private (fld[0:nydim], flux[0:nydim], qIn[0:nydim])
    {

        // advection in y //
        #pragma acc loop gang collapse(2)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                #pragma acc loop vector
                for (int row=0; row<nydim; ++row) {
                    qIn[row] = q2[level][row][col];
                } // end for //
                if (advection_type != 'p') {
                    #pragma acc loop vector
                    for (int row=j1; row<=j2; ++row) {
                        fld[row] = 0.5*(v[level][row][col] + v[level][row+1][col]);
                    } // end for //
                } else {
                    #pragma acc loop vector
                    for (int row=j1; row<=j2+1; ++row) {
                        fld[row] = v[level][row][col];

                        real r = fabs(fld[row]*dt/dy);
                        if (fld[row] >=0.0) {
                            flux[row]= r*(qIn[row-1] + 0.25*(1-r)*(qIn[row]-qIn[row-2]));
                        } else {
                            flux[row]= r*(-qIn[row]  + 0.25*(1-r)*(qIn[row+1] - qIn[row-1]));
                        } // end if //
                    } // end for //
                    #pragma acc loop vector
                    for (int row=j1; row<=j2; ++row) {
                        //qOut[row] = qIn[row] - (flux[row+1] - flux[row]) + dt/dy*qIn[row]*(fld[row+1] - fld[row]);
                        q2[level][row][col] = qIn[row] *(1.0 +  (fld[row+1] - fld[row])*dt/dy ) - (flux[row+1] - flux[row]);
                    } // end for //
                } // end if //
                //advec1(qOut,qIn,fld,flux,dy,dt,j1,j2,nydim,advection_type);
            } // end for //
        } // end for //
    }  // end of parallel region //
    free(qIn);
    free(fld);
    free(flux);


    qIn  = malloc(nzdim*sizeof(real));
    fld  = malloc(nzdim*sizeof(real));
    flux = malloc(nzdim*sizeof(real));

    #pragma acc parallel vector_length(VEC_LEN) present(q2,w,q1) private (fld[0:nzdim], flux[0:nzdim], qIn[0:nzdim])
    {
         // advection in z //
        #pragma acc loop gang collapse(2)
        for (int row=j1; row<=j2; ++row) {
            //#pragma acc loop seq
            for (int col=i1; col<=i2; ++col) {
                #pragma acc loop vector
                for (int level=0; level<nzdim; ++level) {
                    qIn[level] = q2[level][row][col];
                } // end for //
                if (advection_type != 'p') {
                    #pragma acc loop vector
                    for (int level=k1; level<=k2; ++level) {
                        fld[level] = 0.5*(w[level][row][col] + w[level+1][row][col]);
                    } // end for //
                } else {
                    #pragma acc loop vector
                    for (int level=k1; level<=k2+1; ++level) {
                        fld[level] = w[level][row][col];

                        real r = fabs(fld[level]*dt/dz);
                        if (fld[level] >=0.0) {
                            flux[level]= r*(qIn[level-1] + 0.25*(1-r)*(qIn[level]-qIn[level-2]));
                        } else {
                            flux[level]= r*(-qIn[level]  + 0.25*(1-r)*(qIn[level+1] - qIn[level-1]));
                        } // end if //
                    } // end for //
                    #pragma acc loop vector
                    for (int level=k1; level<=k2; ++level) {
                        q2[level][row][col] = qIn[level]*(1.0 + (fld[level+1] - fld[level])*dt/dz ) - (flux[level+1] - flux[level]);
                        //qOut[level] = qIn[level] - (flux[level+1] - flux[level]) + dt/dz*qIn[level]*(fld[level+1] - fld[level]);
                    } // end for //
                } // end if //
                //advec1(qOut,qIn,fld,flux,dz,dt,k1,k2,nzdim,advection_type);
            } // end for //
        } // end for //
    } // end of parallel region //

    free(qIn);
    free(fld);
    free(flux);

    #pragma acc parallel vector_length(VEC_LEN) present(u,v,w,uu,vv,ww)
    {
       	real tempx,tempy,tempz;

        tempx = -tstep/(4*dx);
        tempy = -tstep/(4*dy);
        tempz = -tstep/(4*dz);

        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2; ++row) {
                //#pragma acc loop vector
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

        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2+1; ++row) {
                //#pragma acc loop vector
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

        #pragma acc loop gang collapse(3)
        for (int level=k1; level<=k2+1; ++level) {
            //#pragma acc loop seq
            for (int row=j1; row<=j2; ++row) {
                //#pragma acc loop vector
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
    } // end of parallel region
} // end of advection() //
