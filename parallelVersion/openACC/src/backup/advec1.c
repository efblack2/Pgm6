#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "real.h"

/*
 * ======================= advect1 ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
 *	uv	real array	true speed of wave
 *	dx	real		grid spacing
 *	dt	real		time step
 *	i1,i2	integers	indices bounding array data
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */

#pragma acc routine vector
void advec1(real *restrict q2,const real *restrict q1,real *restrict uv,real *restrict flux,real dx,real dt,int i1,int i2, int nxydim, char advection_type)
{

    real courant;
    int k = 0;
	
    switch ( advection_type ) {
    case 'l':
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i) {
            courant = uv[i] * dt / dx;
            q2[i] = q1[i] - 0.5 * courant * ((q1[i+1] - q1[i-1]) - courant*(q1[i+1] -2.0*q1[i] + q1[i-1]));
        } // end for //
        break;
    case 'c':
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i, ++k) {
            courant = uv[i] * dt / dx;
            q2[i] = q1[i] + 
            courant*                                        (  (q1[i-3]-q1[i+3]) -  9*(q1[i-2]-q1[i+2]) +  45*(q1[i-1]-q1[i+1])            ) /  60.0 +
            courant*courant*                                (2*(q1[i-3]+q1[i+3]) - 27*(q1[i-2]+q1[i+2]) + 270*(q1[i-1]+q1[i+1]) - 490*q1[i]) / 360.0 +
            courant*courant*courant*                        ( -(q1[i-3]-q1[i+3]) +  8*(q1[i-2]-q1[i+2]) -  13*(q1[i-1]-q1[i+1])            ) /  48.0 +
            courant*courant*courant*courant*                ( -(q1[i-3]+q1[i+3]) + 12*(q1[i-2]+q1[i+2]) -  39*(q1[i-1]+q1[i+1]) +  56*q1[i]) / 144.0 +
            courant*courant*courant*courant*courant*        (  (q1[i-3]-q1[i+3]) -  4*(q1[i-2]-q1[i+2]) +   5*(q1[i-1]-q1[i+1])            ) / 240.0 +
            courant*courant*courant*courant*courant*courant*(  (q1[i-3]+q1[i+3]) -  6*(q1[i-2]+q1[i+2]) +  15*(q1[i-1]+q1[i+1]) -  20*q1[i]) / 720.0;
        } // end for //
        break;
    case 't':
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i) {
            courant = uv[i] * dt / dx;
            q2[i] = q1[i] + 0.5*courant*( (q1[i-1] - q1[i+1]) + courant*(q1[i+1] -2.0*q1[i] + q1[i-1]));
            if (courant >= 0) {
                q2[i] -= ( (1+courant)      *courant*(courant-1)/6.0*(q1[i+1] - 3*q1[i] + 3*q1[i-1] - q1[i-2])  );         
            } else {
                q2[i] -= ( (1+fabs(courant))*courant*(courant+1)/6.0*(q1[i-1] - 3*q1[i] + 3*q1[i+1] - q1[i+2])  );    
            } // end if //         
        } // end for //
        break;
    case 'p':
        #pragma acc loop vector
        for (int i=i1; i<=i2+1; ++i) {
            real r = fabs(uv[i]*dt/dx);
            if (uv[i] >=0.0) {
                flux[i]= r*(q1[i-1] + 0.25*(1-r)*(q1[i]-q1[i-2]));
            } else {
                flux[i]= r*(-q1[i]  + 0.25*(1-r)*(q1[i+1] - q1[i-1]));	    
            } // end if //
        } // end for //
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i) {
            q2[i] = q1[i] - (flux[i+1] - flux[i]) + dt/dx*q1[i]*(uv[i+1] - uv[i]);
        } // end for //
        break;
    case 'u':
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i) {
            courant = uv[i] * dt / dx;
            if (courant >= 0.0) {
                q2[i] = q1[i] - courant*(q1[i] - q1[i-1]);
            } else {
                q2[i] = q1[i] - courant*(q1[i+1] - q1[i]);
            } // end if //
        } // end for //
        break;
    case 'n':
        #pragma acc loop vector
        for (int i=i1; i<=i2; ++i,++k) {
            courant = q1[k] * dt / dx;
            q2[i] = q1[i] - 0.5*courant * ((q1[i+1] - q1[i-1]) - courant*(q1[i+1] -2.0*q1[i] + q1[i-1]));
        } // end for //
        break;
    default:
//        printf("Integrate: Error, unrecognized advection type '%c'\n",advection_type);
//        exit(1);
        break;
    } // end switch //

} // end of advec1() //
