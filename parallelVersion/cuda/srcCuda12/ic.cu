#include "real.h"

void getDensity(real *ro_u,real *ro_w,int kB,int kT,real dz,real g,real thetaBar);
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
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		number of grid points
 *
 */

void ic(real ***t ,real ***u ,real ***v ,real ***w ,real *ro_u,real *ro_w,real dx,real dy,real  dz,real uPertur, int iW,int iE,int jS,int jN,int kB,int kT,int bcW,int nx,int ny,int nz,real *x0,real *y0,real *z0, real *deltaTheta,real *deltaV,real thetaBar, real *rx, real *ry, real *rz,real g)
{


    /* initial condition for p field */
	real z=0.5*dz;
    for (int level=kB; level<=kT; ++level, z+=dz) {
        real y=0.5*dy;
        for (int row=jS; row<=jN; ++row, y+=dy) {
            real x=0.5*dx;    
            for (int col=iW; col<=iE; ++col, x+=dx) {
                
                /* initial condition for scalar "V" field */
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        v[level][row][col] = 0.5*deltaV[m] *(cos(rm*M_PI)+1); 
                    } // end if //
                } // end for //
                
                //real const uPertur=DELTAU;
                u[level][row][col] =  uPertur*( (real) rand()/(RAND_MAX + 1.0)) -0.5*uPertur;
                /* uPertur*((real) rand()/(RAND_MAX + 1.0)) - upertur*0.5 ;*/
                /* initial condition for scalar "V" field */
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        v[level][row][col] = 0.5*deltaV[m] *(cos(rm*M_PI)+1); 
                    } // end if //
                } // end for //
            } // end for //
            u[level][row][iE+1] = u[level][row][iE];
        } // end for //
    } // end for //
    
    for (int level=kB; level<=kT; ++level) {
        for (int col=iW; col<=iE; ++col) {
            v[level][jN+1][col] = v[level][jN][col];
        } // end for //
    } // end for //
    
    for (int row=jS; row<=jN; ++row) {
        for (int col=iW; col<=iE; ++col) {
            w[kT+1][row][col] = w[kT][row][col];
        } // end for //
    } // end for //

    
    /* initial condition for scalar "q" field */
	z=0.5*dz;
    for (int level=kB; level<=kT; ++level, z+=dz) {
    	real y=0.5*dy;
        for (int row=jS; row<=jN; ++row, y+=dy) {
        	real x=0.5*dx;
            for (int col=iW; col<=iE; ++col, x+=dx) {
                t[level][row][col] = thetaBar;
                for (int m=0; m<2; m++) {
                    real rm = sqrt( pow( (x-x0[m])/rx[m],2.0) + pow( (y-y0[m])/ry[m],2.0)  + pow((z-z0[m])/rz[m],2.0));
                    if (rm <= 1.0) {
                        t[level][row][col] += 0.5*deltaTheta[m] *(cos(rm*M_PI)+1);
                    } // end if //
                } // end for //
            } // end for //
        } // end for //
    } // end for //
    getDensity(ro_u,ro_w,kB,kT,dz,g,thetaBar);
	return;
} // end of ic() //


