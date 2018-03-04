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
#include <omp.h>
#include "real.h"
#include "dimCube.h"
#include "prototypes.h"


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
    int nthreads;
    double elapsed_time;
    struct timeval tp;
    FILE *ifp;

    real dx,dy,dz,dt;

    /* constants for Computer Problem 5 and 6*/
    const real g = 9.81;
    const real thetaBar = 300.;

	real deltaTheta[2];      /* two temperature perturbations */
	real deltaV[2];          /* two V perturbations */
	real x0[2],y0[2],z0[2];  /* two locations perturbations */
    real rx[2],ry[2],rz[2];  /* two locations perturbations radius */
	real tstep;

// Arrays and other variables //

    real *ro_u, *ro_w;
    real ***t1, ***t2;
    float *tplot;
    real ***p1, ***p2, ***p3;
    real ***u1, ***u2, ***u3;
    real ***v1, ***v2, ***v3;
    real ***w1, ***w2, ***w3;
	real cs,deltau;
	real ku,kv,kw,kt;
    int nx,ny,nz,i1, i2, j1,j2, k1, k2, nxdim,nydim, nzdim, bc_width=0;
    int nstep,nplot;
    char advection_type;

/*
	// setting 1st perturbation case //
    x0[0]= X0_0;
    y0[0]= Y0_0;
    z0[0]= Z0_0;
    deltaTheta[0] = DELTATHETA_0;
	// setting 2nd perturbation case //
    x0[1]= X0_1;
    y0[1]= Y0_1;
    z0[1]= Z0_1;
    deltaTheta[1] = DELTATHETA_1;

    deltaV[0] = DELTAV_0;
    deltaV[1] = DELTAV_1;

    rx[0]= XRADIUS_0;
    rx[1]= XRADIUS_1;
    ry[0]= YRADIUS_0;
    ry[1]= YRADIUS_1;
    rz[0]= ZRADIUS_0;
    rz[1]= ZRADIUS_1;

*/


    /* Parameters and input .................................... */

    if (argc > 1 ) {
        ifp = fopen(argv[1], "r");
        if (ifp == NULL) {
          fprintf(stderr, "Can't open input file %s!\n", argv[1]);
          exit(1);
        }
    } else {
        printf("Use: %s  filename \n", argv[0]);
        printf("Example: %s  data.txt \n", argv[0]);
        printf("advection_scheme can be: L for Lax-Wendroff, C for Crowley, T for Takacs, P for Piecewise, U for Upstream\n");
        exit(0);
    } // endif //

    if (sizeof(real) == 8) {
        printf("Double precision version\n");
    } else {
        printf("Single precision version\n");
    } // end if

    fscanf(ifp, "%d %d %c",&nstep, &nplot, &advection_type);
    fscanf(ifp, "%d %d %d",&nx, &ny, &nz);
    #ifdef DOUBLE
    fscanf(ifp, "%lf %lf %lf %lf",&dx, &dy, &dz, &dt);
    fscanf(ifp, "%lf",&cs);
    fscanf(ifp, "%lf %lf %lf %lf", &ku,&kv,&kw,&kt);
    fscanf(ifp, "%lf",&deltau);
    fscanf(ifp, "%lf %lf %lf",&x0[0], &y0[0], &z0[0]);
    fscanf(ifp, "%lf %lf %lf",&rx[0], &ry[0], &rz[0]);
    fscanf(ifp, "%lf %lf",&deltaTheta[0], &deltaV[0]);
    fscanf(ifp, "%lf %lf %lf",&x0[1], &y0[1], &z0[1]);
    fscanf(ifp, "%lf %lf %lf",&rx[1], &ry[1], &rz[1]);
    fscanf(ifp, "%lf %lf",&deltaTheta[1], &deltaV[1]);
    #else
    fscanf(ifp, "%f %f %f %f",&dx, &dy, &dz, &dt);
    fscanf(ifp, "%f",&cs);
    fscanf(ifp, "%f %f %f %f", &ku,&kv,&kw,&kt);
    fscanf(ifp, "%f",&deltau);
    fscanf(ifp, "%f %f %f",&x0[0], &y0[0], &z0[0]);
    fscanf(ifp, "%f %f %f",&rx[0], &ry[0], &rz[0]);
    fscanf(ifp, "%f %f",&deltaTheta[0], &deltaV[0]);
    fscanf(ifp, "%f %f %f",&x0[1], &y0[1], &z0[1]);
    fscanf(ifp, "%f %f %f",&rx[1], &ry[1], &rz[1]);
    fscanf(ifp, "%f %f",&deltaTheta[1], &deltaV[1]);
    #endif
    fclose(ifp);


    printf("Program 6, ATMS 502/CSE 566, Fall 2016\n");
	printf("Domain (%3d,%3d, %3d)\n",nx,ny,nz);
	printf("dx=%.5f\n", dx);
	printf("dy=%.5f\n", dy);
	printf("dz=%.5f\n", dz);
	printf("dt=%.5f\n", dt);
	printf("cs=%.5f\n", cs);
	printf("KU=%.5f\n", ku);
	printf("KV=%.5f\n", kv);
	printf("KW=%.5f\n", kw);
	printf("KT=%.5f\n", kt);
	printf("deltau=%.5f\n", deltau);

	printf("x0=%.5f\n", x0[0]);
	printf("y0=%.5f\n", y0[0]);
	printf("z0=%.5f\n", z0[0]);

	printf("rx=%.5f\n", rx[0]);
	printf("ry=%.5f\n", ry[0]);
	printf("rz=%.5f\n", rz[0]);

	printf("deltaTheta=%.5f\n", deltaTheta[0]);
	printf("deltaV=%.5f\n", deltaV[0]);


	printf("x1=%.5f\n", x0[1]);
	printf("y1=%.5f\n", y0[1]);
	printf("z1=%.5f\n", z0[1]);

	printf("rx=%.5f\n", rx[1]);
	printf("ry=%.5f\n", ry[1]);
	printf("rz=%.5f\n", rz[1]);

	printf("deltaTheta=%.5f\n", deltaTheta[1]);
	printf("deltaV=%.5f\n", deltaV[1]);


/*
    printf("Domain (%3d,%3d,%3d)\n",NZ,NY,NX);
    printf("dx=%.5f\n", dx);
    printf("dx=%.5f\n", dy);
    printf("dz=%.5f\n", dz);
    printf("dt=%.5f\n", dt);

	printf(" >>> Using left  cone size =  %10.7e %10.7e %10.7e\n",rx[0],ry[0],rz[0]);
	printf(" >>> Using right cone size =  %10.7e %10.7e %10.7e\n",rx[1],ry[1],rz[1]);
*/

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
    i2          = i1+nx-1;
    j1          = bc_width;
    j2          = j1+ny-1;
    k1          = bc_width;
    k2          = k1+nz-1;
    nxdim       = nx+2*bc_width;
    nydim       = ny+2*bc_width;
    nzdim       = nz+2*bc_width;




    t1 = dimCube(nzdim,nydim,nxdim);
    t2 = dimCube(nzdim,nydim,nxdim);
	tplot = malloc( nx * ny * nz * sizeof(float*));

    p1 = dimCube(nz,ny,nx);
    p2 = dimCube(nz,ny,nx);
    p3 = dimCube(nz,ny,nx);


    u1 = dimCube(nzdim,nydim,nxdim+1);
    u2 = dimCube(nzdim,nydim,nxdim+1);
    u3 = dimCube(nzdim,nydim,nxdim+1);

    v1 = dimCube(nzdim,nydim+1,nxdim);
    v2 = dimCube(nzdim,nydim+1,nxdim);
    v3 = dimCube(nzdim,nydim+1,nxdim);


    w1 = dimCube(nzdim+1,nydim,nxdim);
    w2 = dimCube(nzdim+1,nydim,nxdim);
    w3 = dimCube(nzdim+1,nydim,nxdim);


    ro_u=malloc(nzdim*sizeof(real));
    ro_w=malloc((nzdim+1)*sizeof(real));


    // end of allocating memory for arrays and matrices


    gettimeofday(&tp,NULL);
    elapsed_time = -(tp.tv_sec*1.0e6 + tp.tv_usec);

    /*
     * Set and plot the initial condition
     */

    ic(t1,u1,v1,w1,ro_u,ro_w,dx,dy,dz,deltau,i1,i2,j1,j2,k1,k2,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);
    #pragma omp parallel
    {
        bc(u1,v1,w1,i1,i2,j1,j2,k1,k2,bc_width);
        bc4T(t1,t2,i1,i2,j1,j2,k1,k2,bc_width);
        // acctually only the BC of t1 should be copyed to t2
        #pragma omp for
        for (int level=0; level<nzdim; ++level) {
            for (int row=0; row<nydim; ++row) {
                for (int col=0; col<nxdim; ++col) {
                    t2[level][row][col]=t1[level][row][col];
                } // end for //
            } // end for //
        } // end for //

        // Ploting initial conditions //
        #pragma omp for
        for (int col=i1; col<=i2; ++col) {
            for (int row=j1; row<=j2; ++row) {
                for (int level=k1; level<=k2; ++level) {
                    tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (t1[level][row][col] - thetaBar);
                } // end for //
            } // end for //
        } // end for //
        #pragma omp single
        putfield("T",0.0,tplot, nx,ny,nz);

        #pragma omp for
        for (int col=i1; col<=i2; ++col) {
            for (int row=j1; row<=j2; ++row) {
                for (int level=k1; level<=k2; ++level) {
                    tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(u1[level][row][col]+u1[level][row][col+1]));
                } // end for //
            } // end for //
        } // end for //
        #pragma omp single
        putfield("U",0.0,tplot, nx,ny,nz);

        #pragma omp for
        for (int col=i1; col<=i2; ++col) {
            for (int row=j1; row<=j2; ++row) {
                for (int level=k1; level<=k2; ++level) {
                    tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(v1[level][row][col]+v1[level][row+1][col]));
                } // end for //
            } // end for //
        } // end for //
        #pragma omp single
        putfield("V",0.0,tplot, nx,ny,nz);

        #pragma omp for
        for (int col=i1; col<=i2; ++col) {
            for (int row=j1; row<=j2; ++row) {
                for (int level=k1; level<=k2; ++level) {
                    tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(w1[level][row][col]+w1[level+1][row][col]));
                } // end for //
            } // end for //
        } // end for //
        #pragma omp single
        putfield("W",0.0,tplot, nx,ny,nz);

        #pragma omp for
        for (int col=0; col<nx; ++col) {
            for (int row=0; row<ny; ++row) {
                for (int level=0; level<nz; ++level) {
                    tplot[col*nz*ny + row*nz + level] = (float) (p1[level][row][col]);
                } // end for //
            } // end for //
        } // end for //
        #pragma omp single
        putfield("P",0.0,tplot, nx,ny,nz);

        // end of Plotting initial conditions //

        /*
         * .. Integrate .....................................................
         */

        #pragma omp single
        tstep=dt;
        // end of omp single
        for (int n=1 ; n<=nstep; n++) {
            #pragma omp for nowait
	        for (int level=0; level<nzdim; ++level) {
	            for (int row=0; row<nydim; ++row) {
	                for (int col=0; col<=nxdim; ++col) {
	                    u3[level][row][col]=u1[level][row][col];
	                } // end for //
	            } // end for //
            } // end for //

            #pragma omp for nowait
	        for (int level=0; level<nzdim; ++level) {
	            for (int row=0; row<=nydim; ++row) {
	                for (int col=0; col<nxdim; ++col) {
	                    v3[level][row][col]=v1[level][row][col];
	                } // end for //
	            } // end for //
            } // end for //

	        #pragma omp for // using the implicit barrirer to syncromize
	        for (int level=0; level<=nzdim; ++level) {
	            for (int row=0; row<nydim; ++row) {
	                for (int col=0; col<nxdim; ++col) {
	                    w3[level][row][col]=w1[level][row][col];
	                } // end for //
	            } // end for //
            } // end for //

            //  . . . Compute values at next step                           //


            advection(t2,u3,v3,w3,t1,u2,v2,w2,dx,dy,dz,dt,i1,i2,j1,j2,k1,k2,nxdim,nydim,nzdim,tstep,advection_type);

            diffusion(t2,u3,v3,w3,t1,u1,v1,w1,dx*dx,dy*dy,dz*dz,i1,i2,j1,j2,k1,k2,dt,tstep,ku,kv,kw,kt);// KU,KV,KW,KT);

            pgf(p3,u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,i1,i2,j1,j2,k1,k2,tstep, g,thetaBar,cs*cs,bc_width);

            //  . . . Set boundary conditions for T                     //
            bc4T(t2,t1,i1,i2,j1,j2,k1,k2,bc_width);


	        if (n%nplot == 0) {

                #pragma omp for
                for (int col=i1; col<=i2; ++col) {
                    for (int row=j1; row<=j2; ++row) {
                        for (int level=k1; level<=k2; ++level) {
                            tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (t2[level][row][col] - thetaBar);
                        } // end for //
                    } // end for //
                } // end for //
                #pragma omp single
                putfield("T",(float) dt*n,tplot, nx,ny,nz);

                #pragma omp for
                for (int col=i1; col<=i2; ++col) {
                    for (int row=j1; row<=j2; ++row) {
                        for (int level=k1; level<=k2; ++level) {
                            tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(u3[level][row][col]+u3[level][row][col+1]));
                        } // end for //
                    } // end for //
                } // end for //
                #pragma omp single
                putfield("U",(float) dt*n,tplot, nx,ny,nz);

                #pragma omp for
                for (int col=i1; col<=i2; ++col) {
                    for (int row=j1; row<=j2; ++row) {
                        for (int level=k1; level<=k2; ++level) {
                            tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(v3[level][row][col]+v3[level][row+1][col]));
                        } // end for //
                    } // end for //
                } // end for //
                #pragma omp single
                putfield("V",(float) dt*n,tplot, nx,ny,nz);

                #pragma omp for
                for (int col=i1; col<=i2; ++col) {
                    for (int row=j1; row<=j2; ++row) {
                        for (int level=k1; level<=k2; ++level) {
                            tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(w3[level][row][col]+w3[level+1][row][col]));
                        } // end for //
                    } // end for //
                } // end for //
                #pragma omp single
                putfield("W",(float) dt*n,tplot, nx,ny,nz);

                #pragma omp for
                for (int col=0; col<nx; ++col) {
                    for (int row=0; row<ny; ++row) {
                        for (int level=0; level<nz; ++level) {
                            tplot[col*nz*ny + row*nz + level] = (float) (p3[level][row][col]);
                        } // end for //
                    } // end for //
                } // end for //
                #pragma omp single
                putfield("P",(float) dt*n,tplot, nx,ny,nz);
            } // end if //

            //  . . . Do array update at end of time step                   //
            #pragma omp single
            {
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
            } // end of single
        }	// end of time loop n = 1,...,nstep //
        #pragma omp single nowait
        nthreads = omp_get_num_threads();
        // end of omp single
    } // end of parallel region //
    gettimeofday(&tp,NULL);
    elapsed_time += (tp.tv_sec*1.0e6 + tp.tv_usec);
    printf ("\n\nIt tooks %14.6e seconds for %d threads to finish\n", elapsed_time*1.0e-6, nthreads);

    if (sizeof(real) == 8) {
        printf("Double precision version\n");
    } else {
        printf("Single precision version\n");
    } // end if


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

