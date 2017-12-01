#include <iostream>
using namespace std;

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
 *	iW,iE	integers	indices bounding array data
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */

void advec1(real *__restrict__  q2,const real *__restrict__  q1,real *uv,real dx,real dt,int iW,int iE, int nxydim, char advection_type)
{
    real courant;
    //real *flux; 
    
	
    switch ( advection_type ) {
    case 'l':
        for (int i=iW; i<=iE; ++i) {
            courant = uv[i] * dt / dx;
            q2[i] = q1[i] - 0.5 * courant * ((q1[i+1] - q1[i-1]) - courant*(q1[i+1] -2.0*q1[i] + q1[i-1]));
        } // end for //
        break;
    case 'c':
        for (int i=iW, k=0; i<=iE; ++i, ++k) {
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
        for (int i=iW; i<=iE; ++i) {
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
    /*
        {
        real *__restrict__   flux = (real *) malloc( nxydim*sizeof(real) );
        for (int i=iW; i<=iE+1; ++i) {
            real r = fabs(uv[i]*dt/dx);
            if (uv[i] >=0.0) {
                flux[i]= r*(q1[i-1] + 0.25*(1-r)*(q1[i]-q1[i-2]));
            } else {
                flux[i]= r*(-q1[i]  + 0.25*(1-r)*(q1[i+1] - q1[i-1]));	    
            }
        } // end for //
        for (int i=iW; i<=iE; ++i) {
            q2[i] = q1[i] - (flux[i+1] - flux[i]) + dt/dx*q1[i]*(uv[i+1] - uv[i]);
        } // end for 
        free(flux);
        }
        */
        real flux1,flux2,r;
        r = fabs(uv[iW]*dt/dx);

        if (uv[iW] >=0.0) {
            flux1= r*(q1[iW-1] + 0.25*(1-r)*(q1[iW]-q1[iW-2]));
        } else {
            flux1= r*(-q1[iW]  + 0.25*(1-r)*(q1[iW+1] - q1[iW-1]));	    
        } // end if 
        
        for (int i=iW; i<=iE; ++i) {
            r = fabs(uv[i+1]*dt/dx);
            if (uv[i+1] >=0.0) {
                flux2= r*(q1[i] + 0.25*(1-r)*(q1[i+1]-q1[i-1]));
            } else {
                flux2= r*(-q1[i+1]  + 0.25*(1-r)*(q1[i+2] - q1[i]));	    
            } // end if //
            q2[i] = q1[i] - (flux2 - flux1) + dt/dx*q1[i]*(uv[i+1] - uv[i]);
            flux1=flux2;
        } // end for 
        break;
    case 'u':
        for (int i=iW; i<=iE; ++i) {
            courant = uv[i] * dt / dx;
            if (courant >= 0.0) {
                q2[i] = q1[i] - courant*(q1[i] - q1[i-1]);
            } else {
                q2[i] = q1[i] - courant*(q1[i+1] - q1[i]);
            } // end if //
        } // end for //
        break;
    case 'n':
        for (int i=iW, k=0; i<=iE; ++i,++k) {
            courant = q1[k] * dt / dx;
            q2[i] = q1[i] - 0.5*courant * ((q1[i+1] - q1[i-1]) - courant*(q1[i+1] -2.0*q1[i] + q1[i-1]));
        } // end for //
        break;
    default:
        cout << "Integrate: Error, unrecognized advection type " << advection_type << endl;
        exit(-1);
        break;
    } // end switch //

} // end of advec1() //
