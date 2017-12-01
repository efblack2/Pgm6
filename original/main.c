/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm1:  Linear and nonlinear advection
 *  >>>>> PUT YOUR NAME HERE!  Yes, you!! <<<<<
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>

#include "real.h"
#include "dimCube.h"
#include "prototypes.h"
#include "problemParam.h"


/*
 * Definitions
 * NX is the number of physical grid points - Not including 'ghost points'
 * bc_width is the number of ghost points - on each side of the physical grid
 * i1 is the first C index to use - the left-most physical point
 * i2 is the last C index to use - the right-most physical point
 * nxdim is the total array size, including grid points; use for declarations
 * MAXSTEP is the maximum number of time steps to take
 * "c" is [constant] flow speed (e.g. meters/sec) - used only in linear case
 * "dx" is grid spacing: distance (e.g. meters) between grid points in X, Y
 * "name" is a variable containing your name - used to label the plots.
 */



int main(int argc, char *argv[])
{

    double elapsed_time;
    struct timeval tp;

    real dx   = DX;
    real dy   = DY;
    real dz   = DZ;
    //char *name = "Edgar F. Black";
    const real dt = DT;

    /* constants for Computer Problem 5 and 6*/
    const real g = 9.81;
    const real thetaBar = 300.;
	real deltau=0;
    
    
	real deltaTheta[2];      /* two temperature perturbations */
	real deltaV[2];          /* two V perturbations */
	real x0[2],y0[2],z0[2];  /* two locations perturbations */
    real rx[2],ry[2],rz[2];  /* two locations perturbations radius */


	real tstep;
    
	/* setting 1st perturbation case */
    x0[0]= X0_0;
    y0[0]= Y0_0;
    z0[0]= Z0_0;
    deltaTheta[0] = DELTATHETA_0;
	/* setting 2nd perturbation case */
    x0[1]= X0_1;
    y0[1]= Y0_1;
    z0[1]= Z0_1;
    deltaTheta[1] = DELTATHETA_1;

    #if (defined DELTAV_0 && defined DELTAV_1)
    deltaV[0] = DELTAV_0;
    deltaV[1] = DELTAV_1;
    #endif

    rx[0]= XRADIUS_0;
    rx[1]= XRADIUS_1;
    ry[0]= YRADIUS_0;
    ry[1]= YRADIUS_1;
    rz[0]= ZRADIUS_0;
    rz[1]= ZRADIUS_1;


    
    const real cs2=CS*CS;

    
    /* Arrays and other variables */
    int i1, i2, j1,j2, k1, k2, nxdim,nydim, nzdim, bc_width=0;
    int nstep,nplot;
    
    char advection_type;



    /* Parameters and input .................................... */

    if (argc == 4 ) {
        nstep = atoi(argv[1]);
        nplot = atoi(argv[2]);
        advection_type= tolower(argv[3][0])  ;
    } else {
        printf("Use: %s  'nstep' 'nplot' 'advection_scheme' \n", argv[0]);       
        printf("Example: %s  100 10 l \n", argv[0]);       
        printf("advection_scheme can be: L for Lax-Wendroff, C for Crowley, T for Takacs, P for Piecewise, U for Upstream\n"); exit(0);        
    } // endif //


    printf("Program 6, ATMS 502/CSE 566, Fall 2016\n");
    printf("Domain (%3d,%3d,%3d)\n",NZ,NY,NX);
    printf("dx=%.5f\n", dx);
    printf("dx=%.5f\n", dy);
    printf("dz=%.5f\n", dz);
    printf("dt=%.5f\n", dt);
	
	printf(" >>> Using left  cone size =  %10.7e %10.7e %10.7e\n",rx[0],ry[0],rz[0]);
	printf(" >>> Using right cone size =  %10.7e %10.7e %10.7e\n",rx[1],ry[1],rz[1]);
	
	
    switch ( advection_type ) {
    case 'l':
    case 'u':
        bc_width = 1;
        break;
    case 'p':
    case 't':
        bc_width = 2;
        break;
    case 'c':
        bc_width = 3;
        break;
    case 'n':
        break;
    default:
        printf("Integrate: Error, unrecognized advection type '%c'\n",advection_type);
        exit(1);
        break;
    } // end switch //
    i1          = bc_width;
    i2          = i1+NX-1;
    j1          = bc_width;
    j2          = j1+NY-1;
    k1          = bc_width;
    k2          = k1+NZ-1;
    nxdim       = NX+2*bc_width;
    nydim       = NY+2*bc_width;
    nzdim       = NZ+2*bc_width;


// Arrays and other variables //

    real ***t1, ***t2;
    float *tplot;
    real ***p1, ***p2, ***p3;
    real ***u1, ***u2, ***u3;
    real ***v1, ***v2, ***v3;
    real ***w1, ***w2, ***w3;
    
    
    t1 = dimCube(nzdim,nydim,nxdim);
    t2 = dimCube(nzdim,nydim,nxdim);
	tplot = malloc( NX * NY * NZ * sizeof(float*));

    p1 = dimCube(NZ,NY,NX);
    p2 = dimCube(NZ,NY,NX);
    p3 = dimCube(NZ,NY,NX);


    u1 = dimCube(nzdim,nydim,nxdim+1);
    u2 = dimCube(nzdim,nydim,nxdim+1);
    u3 = dimCube(nzdim,nydim,nxdim+1);
    
    v1 = dimCube(nzdim,nydim+1,nxdim);
    v2 = dimCube(nzdim,nydim+1,nxdim);
    v3 = dimCube(nzdim,nydim+1,nxdim);
    

    w1 = dimCube(nzdim+1,nydim,nxdim);
    w2 = dimCube(nzdim+1,nydim,nxdim);
    w3 = dimCube(nzdim+1,nydim,nxdim);
    

    real *ro_u, *ro_w;
    ro_u=malloc(nzdim*sizeof(real));
    ro_w=malloc((nzdim+1)*sizeof(real));
    
    
    // end of allocating memory for arrays and matrices


    gettimeofday(&tp,NULL);
    elapsed_time = -(tp.tv_sec*1.0e6 + tp.tv_usec);  

    /*
     * Set and plot the initial condition
     */
    
    
    #if (defined TEST_F || defined OFFICIAL) 
    deltau=DELTAU;
    #endif
    ic(t1,p1,u1,v1,w1,ro_u,ro_w,dx,dy,dz,deltau,i1,i2,j1,j2,k1,k2,bc_width,NX,NY,NZ,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);
    bc(u1,v1,w1,i1,i2,j1,j2,k1,k2,bc_width);
    bc4T(t1,i1,i2,j1,j2,k1,k2,bc_width);
    
    // acctually only the BC of t1 should be copyed to t2 
    for (int level=0; level<nzdim; ++level) {
        for (int row=0; row<nydim; ++row) {
            for (int col=0; col<nxdim; ++col) {
                t2[level][row][col]=t1[level][row][col];
            } // end for //
        } // end for //
    } // end for //
    
    
    // Ploting initial conditions //
    for (int col=i1,l=0; col<=i2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[l++] = (float) (t1[level][row][col] - thetaBar);
            } // end for //
        } // end for //
    } // end for // 

    putfield("T",0.0,tplot, NX,NY,NZ);
    
    for (int col=i1,l=0; col<=i2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[l++] = (float) (0.5*(u1[level][row][col]+u1[level][row][col+1]));
            } // end for //
        } // end for //
    } // end for // 
    putfield("U",0.0,tplot, NX,NY,NZ);

    for (int col=i1,l=0; col<=i2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[l++] = (float) (0.5*(v1[level][row][col]+v1[level][row+1][col]));
            } // end for //
        } // end for //
    } // end for //     
    putfield("V",0.0,tplot, NX,NY,NZ);
    
    for (int col=i1,l=0; col<=i2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[l++] = (float) (0.5*(w1[level][row][col]+w1[level+1][row][col]));
            } // end for //
        } // end for //
    } // end for //
    putfield("W",0.0,tplot, NX,NY,NZ);

    for (int col=0,l=0; col<NX; ++col) {
        for (int row=0; row<NY; ++row) {
            for (int level=0; level<NZ; ++level) {
                tplot[l++] = (float) (p1[level][row][col]);
            } // end for //
        } // end for //
    } // end for //
    putfield("P",0.0,tplot, NX,NY,NZ);

    // end of Plotting initial conditions //


    /*
     * .. Integrate .....................................................
     */
    
     
    tstep=dt;     
    for (int n=1 ; n<=nstep; n++) {
    
	    for (int level=0; level<nzdim; ++level) {
	        for (int row=0; row<nydim; ++row) {
	            for (int col=0; col<=nxdim; ++col) {
	                u3[level][row][col]=u1[level][row][col];
	            } // end for //
	        } // end for //
        } // end for //


	    for (int level=0; level<nzdim; ++level) {
	        for (int row=0; row<=nydim; ++row) {
	            for (int col=0; col<nxdim; ++col) {
	                v3[level][row][col]=v1[level][row][col];
	            } // end for //
	        } // end for //
        } // end for //
	    
	    for (int level=0; level<=nzdim; ++level) {
	        for (int row=0; row<nydim; ++row) {
	            for (int col=0; col<nxdim; ++col) {
	                w3[level][row][col]=w1[level][row][col];
	            } // end for //
	        } // end for //
        } // end for //
        
        #if (defined TEST_B || defined TEST_C)
	    for (int level=0; level<nzdim; ++level) {
	        for (int row=0; row<nydim; ++row) {
	            for (int col=0; col<nxdim; ++col) {
	                t2[level][row][col]=t1[level][row][col];
	            } // end for //
	        } // end for //
        } // end for //
        #endif
        //  . . . Compute values at next step                           //

        #if (defined TEST_A1 || defined TEST_A2 || defined TEST_A3  || defined TEST_D1  || defined TEST_D2 || defined TEST_E|| defined TEST_F || defined OFFICIAL) 
        advection(t2,u3,v3,w3,t1,u2,v2,w2,dx,dy,dz,dt,i1,i2,j1,j2,k1,k2,nxdim,nydim,nzdim,tstep,advection_type);
	    #endif

        #if ( defined TEST_C ||  defined TEST_D1 || defined TEST_D2 || defined TEST_E || defined TEST_F || defined OFFICIAL)
        diffusion(t2,u3,v3,w3,t1,u1,v1,w1,dx*dx,dy*dy,dz*dz,i1,i2,j1,j2,k1,k2,dt,tstep,KU,KV,KW,KT);
        #endif

        #if (defined TEST_B || defined TEST_D1 || defined TEST_D2 || defined TEST_E|| defined TEST_F ||defined OFFICIAL)
        pgf(p3,u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,i1,i2,j1,j2,k1,k2,tstep, g,thetaBar,cs2,bc_width);
	    #endif

        //  . . . Set boundary conditions for T                     //
        bc4T(t2,i1,i2,j1,j2,k1,k2,bc_width);


	    if (n%nplot == 0) {
	    

            for (int col=i1,l=0; col<=i2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[l++] = (float) (t2[level][row][col] - thetaBar);
                    } // end for //
                } // end for //
            } // end for // 
            putfield("T",(float) dt*n,tplot, NX,NY,NZ);
            
            for (int col=i1,l=0; col<=i2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[l++] = (float) (0.5*(u3[level][row][col]+u3[level][row][col+1]));
                    } // end for //
                } // end for //
            } // end for // 
            putfield("U",(float) dt*n,tplot, NX,NY,NZ);

            for (int col=i1,l=0; col<=i2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[l++] = (float) (0.5*(v3[level][row][col]+v3[level][row+1][col]));
                    } // end for //
                } // end for //
            } // end for //     
            putfield("V",(float) dt*n,tplot, NX,NY,NZ);
            
            for (int col=i1,l=0; col<=i2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[l++] = (float) (0.5*(w3[level][row][col]+w3[level+1][row][col]));
                    } // end for //
                } // end for //
            } // end for //
            putfield("W",(float) dt*n,tplot, NX,NY,NZ);

            for (int col=0,l=0; col<NX; ++col) {
                for (int row=0; row<NY; ++row) {
                    for (int level=0; level<NZ; ++level) {
                        tplot[l++] = (float) (p3[level][row][col]);
                    } // end for //
                } // end for //
            } // end for //
            putfield("P",(float) dt*n,tplot, NX,NY,NZ);
	    
        } // end if //

        //  . . . Do array update at end of time step                   //
        if (n == 1) { 
            tstep=2.0*dt;
            update(&p2,&p3);
            update(&u2,&u3);
            update(&v2,&v3);
            update(&w2,&w3);
            
        } else {
            update(&p1,&p2);
            update(&p2,&p3);
            
            update(&u1,&u2);
            update(&u2,&u3);
            
            update(&v1,&v2);
            update(&v2,&v3);
            
            update(&w1,&w2);
            update(&w2,&w3);
        } // end if //

        update(&t1,&t2);

    }	// end of time loop n = 1,...,nstep //
    gettimeofday(&tp,NULL);
    elapsed_time += (tp.tv_sec*1.0e6 + tp.tv_usec);
    printf ("\n\nIt tooks %14.6e seconds to finish\n", elapsed_time*1.0e-6);

    
    freeCube(&t2);
    freeCube(&t1);

    freeCube(&p1);
    freeCube(&p2);
    freeCube(&p3);
    
    freeCube(&u1);
    freeCube(&u2);
    freeCube(&u3);

	freeCube(&v1);
	freeCube(&v2);
	freeCube(&v3);

    freeCube(&w1);
    freeCube(&w2);
    freeCube(&w3);

	free(tplot);
	free(ro_u);
	free(ro_w);

    
    return 0;
} // end main() //

