#include <stdlib.h>
#include <math.h>
#include "real.h"
#include "problemParam.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

void getDensity(real *ro_u,real *ro_w,int k1,int k2,real dz,real g,real thetaBar);
/*
 * ============================ ic =====================
 * IC sets the initial condition
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	IC data. Set 1..nx here;
 *				  [0],[nx+1] = ghost zones
 *				  if 1 ghost point on each side
 *	dx	real		grid spacing
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of grid points
 *
 */

void ic(real ***t ,real ***p ,real ***u ,real ***v ,real ***w ,real *ro_u,real *ro_w,real dx,real dy,real  dz,real uPertur, int i1,int i2,int j1,int j2,int k1,int k2,int bcW,int nx,int ny,int nz,real *x0,real *y0,real *z0, real *deltaTheta,real *deltaV,real thetaBar, real *rx, real *ry, real *rz,real g)
{


    /* initial condition for p field */
	real z=0.5*dz;
    for (int level=k1; level<=k2; ++level, z+=dz) {
        real y=0.5*dy;
        for (int row=j1; row<=j2; ++row, y+=dy) {
            real x=0.5*dx;    
            for (int col=i1; col<=i2; ++col, x+=dx) {
                #if defined TEST_A1
                u[level][row][col] = 20.0;
                v[level][row][col] = 0.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                #elif defined TEST_A2
                u[level][row][col] = 0.0;
                v[level][row][col] = 20.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                #elif defined TEST_A3
                u[level][row][col] = 0.0;
                v[level][row][col] = 0.0;
                w[level][row][col] = -20.0;
                //p[level][row][col] = 0.0;
                #elif defined TEST_B
                u[level][row][col] = 0.0;
                v[level][row][col] = 0.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                #elif defined TEST_C
                u[level][row][col] = 0.0;
                v[level][row][col] = 0.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                #elif (defined TEST_D1 || defined TEST_D2 )
                u[level][row][col] = 0.0;
                v[level][row][col] = 0.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                #elif (defined TEST_E )
                u[level][row][col] = 0.0;
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                /* initial condition for scalar "V" field */
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        v[level][row][col] = 0.5*deltaV[m] *(cos(rm*M_PI)+1); 
                        /*printf("i,j,k, deltaV: %d %d %d %f\n", i,j,k,deltaV[m]);*/
                    } // end if //
                } // end for //
                #elif (defined TEST_F || defined OFFICIAL)
                //real const uPertur=DELTAU;
                u[level][row][col] =  uPertur*( (real) rand()/(RAND_MAX + 1.0)) -0.5*uPertur;
                /* uPertur*((real) rand()/(RAND_MAX + 1.0)) - upertur*0.5 ;*/
                w[level][row][col] = 0.0;
                //p[level][row][col] = 0.0;
                /* initial condition for scalar "V" field */
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        v[level][row][col] = 0.5*deltaV[m] *(cos(rm*M_PI)+1); 
                        /*printf("i,j,k, deltaV: %d %d %d %f\n", i,j,k,deltaV[m]);*/
                    } // end if //
                } // end for //
                #endif
            } // end for //
            u[level][row][i2+1] = u[level][row][i2];
        } // end for //
    } // end for //
    
    for (int level=k1; level<=k2; ++level) {
        for (int col=i1; col<=i2; ++col) {
            v[level][j2+1][col] = v[level][j2][col];
        } // end for //
    } // end for //
    
    for (int row=j1; row<=j2; ++row) {
        for (int col=i1; col<=i2; ++col) {
            w[k2+1][row][col] = w[k2][row][col];
        } // end for //
    } // end for //
    
/*
    for (i=i1; i<=i2; i++) {
        w[k2+1][col] = w[k2][col];
    }    
*/

    
    /* initial condition for scalar "q" field */
	z=0.5*dz;
    for (int level=k1; level<=k2; ++level, z+=dz) {
    	real y=0.5*dy;
        for (int row=j1; row<=j2; ++row, y+=dy) {
        	real x=0.5*dx;
            for (int col=i1; col<=i2; ++col, x+=dx) {
                t[level][row][col] = thetaBar;
                //printf("t[][][]: %f\n", t[level][row][col]);
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        t[level][row][col] += 0.5*deltaTheta[m] *(cos(rm*M_PI)+1);
                        //printf("m, deltatheta: %d, %f %f %f %f\n", m, deltaTheta[m], x, z,t[level][row][col]); 
                    } // end if //
                } // end for //
            } // end for //
        } // end for //
    } // end for //
    getDensity(ro_u,ro_w,k1,k2,dz,g,thetaBar);
	return;
} // end of ic() //


