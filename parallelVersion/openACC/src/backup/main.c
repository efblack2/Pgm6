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
    real ***restrict t1, ***restrict t2;
    float *tplot;
    real ***restrict p1, ***restrict p2, ***restrict p3;
    real ***restrict u1, ***restrict u2, ***restrict u3;
    real ***restrict v1, ***restrict v2, ***restrict v3;
    real ***restrict w1, ***restrict w2, ***restrict w3;
	real cs,deltau;
	real ku,kv,kw,kt;
    int nx,ny,nz,i1, i2, j1,j2, k1, k2, nxdim,nydim, nzdim, bc_width=0;
    int nstep,nplot;
    char advection_type;
    int myVoid;

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

    myVoid=fscanf(ifp, "%d %d %c",&nstep, &nplot, &advection_type);
    myVoid=fscanf(ifp, "%d %d %d",&nx, &ny, &nz);
    #ifdef FLOAT     
    myVoid=fscanf(ifp, "%f %f %f %f",&dx, &dy, &dz, &dt);
    myVoid=fscanf(ifp, "%f",&cs);
    myVoid=fscanf(ifp, "%f %f %f %f", &ku,&kv,&kw,&kt);
    myVoid=fscanf(ifp, "%f",&deltau);
    myVoid=fscanf(ifp, "%f %f %f",&x0[0], &y0[0], &z0[0]);
    myVoid=fscanf(ifp, "%f %f %f",&rx[0], &ry[0], &rz[0]);
    myVoid=fscanf(ifp, "%f %f",&deltaTheta[0], &deltaV[0]);
    myVoid=fscanf(ifp, "%f %f %f",&x0[1], &y0[1], &z0[1]);
    myVoid=fscanf(ifp, "%f %f %f",&rx[1], &ry[1], &rz[1]);
    myVoid=fscanf(ifp, "%f %f",&deltaTheta[1], &deltaV[1]);
    #else
    myVoid=fscanf(ifp, "%lf %lf %lf %lf",&dx, &dy, &dz, &dt);
    myVoid=fscanf(ifp, "%lf",&cs);
    myVoid=fscanf(ifp, "%lf %lf %lf %lf", &ku,&kv,&kw,&kt);
    myVoid=fscanf(ifp, "%lf",&deltau);
    myVoid=fscanf(ifp, "%lf %lf %lf",&x0[0], &y0[0], &z0[0]);
    myVoid=fscanf(ifp, "%lf %lf %lf",&rx[0], &ry[0], &rz[0]);
    myVoid=fscanf(ifp, "%lf %lf",&deltaTheta[0], &deltaV[0]);
    myVoid=fscanf(ifp, "%lf %lf %lf",&x0[1], &y0[1], &z0[1]);
    myVoid=fscanf(ifp, "%lf %lf %lf",&rx[1], &ry[1], &rz[1]);
    myVoid=fscanf(ifp, "%lf %lf",&deltaTheta[1], &deltaV[1]);
    #endif



    printf("Program 6, ATMS 502/CSE 566, Fall 2016\n");
	printf("Domain (%3d,%3d, %3d)\n",nx,ny,nz);
	printf("dx=%.5f\n", dx);
	printf("dy=%.5f\n", dy);
	printf("dz=%.5f\n", dz);
	printf("dt=%.5f\n", dt);
	printf("cs=%.5f\n", cs);
	printf("KU=%.5f\n", ku);
	printf("KV=%.5f\n", kv);
	printf("KV=%.5f\n", kw);
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
     
    ic(t1,p1,u1,v1,w1,ro_u,ro_w,dx,dy,dz,deltau,i1,i2,j1,j2,k1,k2,bc_width,nx,ny,nz,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);

   
    #pragma acc data    pcopyin(t1[0:nzdim][0:nydim][0:nxdim]),\
                        pcreate(t2[0:nzdim][0:nydim][0:nxdim]),\
                        pcopyin(p1[0:nz][0:ny][0:nx]),\
                        pcopyin(p2[0:nz][0:ny][0:nx]),\
                        pcopyin(p3[0:nz][0:ny][0:nx]),\
                        pcopyin(u1[0:nzdim][0:nydim][0:nxdim+1]),\
                        pcopyin(u2[0:nzdim][0:nydim][0:nxdim+1]),\
                        pcopyin(u3[0:nzdim][0:nydim][0:nxdim+1]),\
                        pcopyin(v1[0:nzdim][0:nydim+1][0:nxdim]),\
                        pcopyin(v2[0:nzdim][0:nydim+1][0:nxdim]),\
                        pcopyin(v3[0:nzdim][0:nydim+1][0:nxdim]),\
                        pcopyin(w1[0:nzdim+1][0:nydim][0:nxdim]),\
                        pcopyin(w2[0:nzdim+1][0:nydim][0:nxdim]),\
                        pcopyin(w3[0:nzdim+1][0:nydim][0:nxdim]),\
                        pcopyin(ro_u[0:nzdim]),\
                        pcopyin(ro_w[0:nzdim+1])
    {
    
        bc(u1,v1,w1,i1,i2,j1,j2,k1,k2,bc_width);
        bc4T(t1,t2,i1,i2,j1,j2,k1,k2,bc_width);
        
        
        // Ploting initial conditions //
        //#pragma acc update self(t1[0:nzdim][0:nydim][0:nxdim],\
                                u1[0:nzdim][0:nydim][0:nxdim+1],\
                                v1[0:nzdim][0:nydim+1][0:nxdim],\
                                w1[0:nzdim+1][0:nydim][0:nxdim])

        #pragma omp parallel
        {
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
        } // end of omp parallel region
        // end of Plotting initial conditions //
         
        /*
         * .. Integrate .....................................................
         */
        
        tstep=dt;     
        for (int n=1; n<=nstep; ++n) {
            
            #pragma acc parallel present(u3[0:nzdim][0:nydim][0:nxdim+1],\
                                         u1[0:nzdim][0:nydim][0:nxdim+1],\
                                         v3[0:nzdim][0:nydim+1][0:nxdim],\
                                         v1[0:nzdim][0:nydim+1][0:nxdim],\
                                         w3[0:nzdim+1][0:nydim][0:nxdim],\
                                         w1[0:nzdim+1][0:nydim][0:nxdim])

//                                       t2[0:nzdim][0:nydim][0:nxdim],\
                                         t1[0:nzdim][0:nydim][0:nxdim],\

                                         
            {
                #pragma acc loop gang collapse(3)
                for (int level=0; level<nzdim; ++level) {
                    //#pragma acc loop seq
                    for (int row=0; row<nydim; ++row) {
                        //#pragma acc loop vector 
                        for (int col=0; col<=nxdim; ++col) {
                            u3[level][row][col]=u1[level][row][col];
                        } // end for //
                    } // end for //
                } // end for //

                #pragma acc loop gang collapse(3)
                for (int level=0; level<nzdim; ++level) {
                    //#pragma acc loop seq
                    for (int row=0; row<=nydim; ++row) {
                        //#pragma acc loop vector 
                        for (int col=0; col<nxdim; ++col) {
                            v3[level][row][col]=v1[level][row][col];
                        } // end for //
                    } // end for //
                } // end for //

                #pragma acc loop gang collapse(3)
                for (int level=0; level<=nzdim; ++level) {
                    //#pragma acc loop seq
                    for (int row=0; row<nydim; ++row) {
                        //#pragma acc loop vector 
                        for (int col=0; col<nxdim; ++col) {
                            w3[level][row][col]=w1[level][row][col];
                        } // end for //
                    } // end for //
                } // end for //              
            } // end of parallel region

            
            //  . . . Compute values at next step                           //

            advection(t2,u3,v3,w3,t1,u2,v2,w2,dx,dy,dz,dt,i1,i2,j1,j2,k1,k2,nxdim,nydim,nzdim,tstep,advection_type);

            diffusion(t2,u3,v3,w3,t1,u1,v1,w1,dx*dx,dy*dy,dz*dz,i1,i2,j1,j2,k1,k2,dt,tstep,ku,kv,kw,kt);// KU,KV,KW,KT);
            
            pgf(p3,u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,i1,i2,j1,j2,k1,k2,tstep, g,thetaBar,cs*cs,bc_width);

            //  . . . Set boundary conditions for T                     //
            bc4T(t2,t1,i1,i2,j1,j2,k1,k2,bc_width);

	        if (n%nplot == 0) {
	            #pragma acc update self(t2[0:nzdim][0:nydim][0:nxdim],\
	                                    u3[0:nzdim][0:nydim][0:nxdim+1],\
	                                    v3[0:nzdim][0:nydim+1][0:nxdim],\
	                                    w3[0:nzdim+1][0:nydim][0:nxdim],\
	                                    p3[0:nz][0:ny][0:nx]) 
	         
	            #pragma omp parallel
	            {
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
	            } // end of parallel omp region //
            } // end if //

            //  . . . Do array update at end of time step                   //
            if (n == 1) { 
                tstep=2.0*dt;
                
                #pragma acc parallel vector_length(128)  present(\
                                                         p2[0:nz][0:ny][0:nx],\
                                                         p3[0:nz][0:ny][0:nx],\
                                                         u2[0:nzdim][0:nydim][0:nxdim+1],\
                                                         u3[0:nzdim][0:nydim][0:nxdim+1],\
                                                         v2[0:nzdim][0:nydim+1][0:nxdim],\
                                                         v3[0:nzdim][0:nydim+1][0:nxdim],\
                                                         w2[0:nzdim+1][0:nydim][0:nxdim],\
                                                         w3[0:nzdim+1][0:nydim][0:nxdim],\
                                                         t1[0:nzdim][0:nydim][0:nxdim],\
                                                         t2[0:nzdim][0:nydim][0:nxdim])
                {                
	                real **mtx;
                    #pragma acc loop gang, vector
                    for(int k=0; k < nz; ++k){
                        mtx = p2[k];
                        p2[k] = p3[k];
                        p3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = u2[k];
                        u2[k] = u3[k];
                        u3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = v2[k];
                        v2[k] = v3[k];
                        v3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k <=nzdim; ++k){
                        mtx = w2[k];
                        w2[k] = w3[k];
                        w3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = t2[k];
                        t2[k] = t1[k];
                        t1[k] = mtx;
                    } // end for //
                } // end of parallel region //                
            } else {            
                #pragma acc parallel vector_length(128)  present(\
                                                         p1[0:nz][0:ny][0:nx],\
                                                         p2[0:nz][0:ny][0:nx],\
                                                         p3[0:nz][0:ny][0:nx],\
                                                         u1[0:nzdim][0:nydim][0:nxdim+1],\
                                                         u2[0:nzdim][0:nydim][0:nxdim+1],\
                                                         u3[0:nzdim][0:nydim][0:nxdim+1],\
                                                         v1[0:nzdim][0:nydim+1][0:nxdim],\
                                                         v2[0:nzdim][0:nydim+1][0:nxdim],\
                                                         v3[0:nzdim][0:nydim+1][0:nxdim],\
                                                         w1[0:nzdim+1][0:nydim][0:nxdim],\
                                                         w2[0:nzdim+1][0:nydim][0:nxdim],\
                                                         w3[0:nzdim+1][0:nydim][0:nxdim],\
                                                         t1[0:nzdim][0:nydim][0:nxdim],\
                                                         t2[0:nzdim][0:nydim][0:nxdim])
                                                         
                
                {                
	                real **mtx;
                    #pragma acc loop gang, vector
                    for(int k=0; k < nz; ++k){
                        mtx = p2[k];
                        p2[k] = p1[k];
                        p1[k] = mtx;
                    } // end for //
                    #pragma acc loop gang, vector
                    for(int k=0; k < nz; ++k){
                        mtx = p2[k];
                        p2[k] = p3[k];
                        p3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = u2[k];
                        u2[k] = u1[k];
                        u1[k] = mtx;
                    } // end for //
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = u2[k];
                        u2[k] = u3[k];
                        u3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = v2[k];
                        v2[k] = v1[k];
                        v1[k] = mtx;
                    } // end for //
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = v2[k];
                        v2[k] = v3[k];
                        v3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k <= nzdim; ++k){
                        mtx = w2[k];
                        w2[k] = w1[k];
                        w1[k] = mtx;
                    } // end for //
                    #pragma acc loop gang, vector
                    for(int k=0; k <= nzdim; ++k){
                        mtx = w2[k];
                        w2[k] = w3[k];
                        w3[k] = mtx;
                    } // end for //
                    
                    #pragma acc loop gang, vector
                    for(int k=0; k < nzdim; ++k){
                        mtx = t2[k];
                        t2[k] = t1[k];
                        t1[k] = mtx;
                    } // end for //
                } // end of parallel region //                 
            } // end if //
            
        }	// end of time loop n = 1,...,nstep //
    } // end of data region //
    gettimeofday(&tp,NULL);
    elapsed_time += (tp.tv_sec*1.0e6 + tp.tv_usec);
    printf ("\n\nIt tooks %14.6e seconds to finish\n", elapsed_time*1.0e-6);

    
    freeCube((real ****) &t2);
    freeCube((real ****) &t1);

    freeCube((real ****) &p1);
    freeCube((real ****) &p2);
    freeCube((real ****) &p3);
    
    freeCube((real ****) &u1);
    freeCube((real ****) &u2);
    freeCube((real ****) &u3);

	freeCube((real ****) &v1);
	freeCube((real ****) &v2);
	freeCube((real ****) &v3);

    freeCube((real ****) &w1);
    freeCube((real ****) &w2);
    freeCube((real ****) &w3);

	free(tplot);
	free(ro_u);
	free(ro_w);

    
    return 0;
} // end main() //

