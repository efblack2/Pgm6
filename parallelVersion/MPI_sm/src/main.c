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
#include <mpi.h>
#include "real.h"
#include "dimCube.h"
#include "prototypes.h"
#include "myMPI.h"

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

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    int myRank, commSize;

    MPI_Comm sm_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED, 0,MPI_INFO_NULL, &sm_comm);
    MPI_Comm_rank(sm_comm,&myRank);
    MPI_Comm_size(sm_comm,&commSize);
    const int root=0;
    //const int root=commSize-1;

    //int nthreads;
    double elapsed_time;

    FILE *ifp=NULL;

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

    /* Parameters and input .................................... */

    if (argc > 1 ) {
        ifp = fopen(argv[1], "r");
        if (ifp == NULL) {
            if (myRank == root) fprintf(stderr, "Can't open input file %s!\n", argv[1]);
            MPI_Finalize();
            exit(1);
        } // end if //
    } else {
        if (myRank == root) {
            printf("Use: %s  filename \n", argv[0]);
            printf("Example: %s  data.txt \n", argv[0]);
            printf("advection_scheme can be: L for Lax-Wendroff, C for Crowley, T for Takacs, P for Piecewise, U for Upstream\n");
        } // end if //
        MPI_Finalize();
        exit(0);
    } // endif //

    if (myRank == root) {
        if (sizeof(real) == 8) {
            printf("Double precision version\n");
        } else {
            printf("Single precision version\n");
        } // end if
    } // end if

    int file_free = 0;
    MPI_Status status;

    if (myRank == 0) {
        file_free = 1;
    } else {
        MPI_Recv(&file_free, 1, MPI_INT,myRank-1,1,sm_comm,&status);
    } // end if //

    if (file_free == 1) {
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
    } // end if //
    fclose(ifp);

    if (myRank != commSize-1) {
        MPI_Send(&file_free, 1, MPI_INT,myRank+1,1,sm_comm);
    } // end if //



    if (myRank == root) {
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
    } // end if //

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
        if (myRank == root) printf("Integrate: Error, unrecognized advection type '%c'\n",advection_type);
        MPI_Finalize();
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

/*
    MPI_Aint sz;
    int dispUnit;
*/


    MPI_Win sm_win_t1;
    MPI_Win sm_win_t2;

    t1 = dimCube(nzdim,nydim,nxdim,bc_width, &sm_win_t1, &sm_comm);
    t2 = dimCube(nzdim,nydim,nxdim,bc_width, &sm_win_t2, &sm_comm);

    MPI_Win sm_win_p1;
    MPI_Win sm_win_p2;
    MPI_Win sm_win_p3;

    p1 = dimCube(nz,ny,nx,0, &sm_win_p1, &sm_comm);
    p2 = dimCube(nz,ny,nx,0, &sm_win_p2, &sm_comm);
    p3 = dimCube(nz,ny,nx,0, &sm_win_p3, &sm_comm);

    MPI_Win sm_win_u1;
    MPI_Win sm_win_u2;
    MPI_Win sm_win_u3;

    u1 = dimCube(nzdim,nydim,(nxdim+1),bc_width, &sm_win_u1, &sm_comm);
    u2 = dimCube(nzdim,nydim,(nxdim+1),bc_width, &sm_win_u2, &sm_comm);
    u3 = dimCube(nzdim,nydim,(nxdim+1),bc_width, &sm_win_u3, &sm_comm);

    MPI_Win sm_win_v1;
    MPI_Win sm_win_v2;
    MPI_Win sm_win_v3;

    v1 = dimCube(nzdim,(nydim+1),nxdim,bc_width, &sm_win_v1, &sm_comm);
    v2 = dimCube(nzdim,(nydim+1),nxdim,bc_width, &sm_win_v2, &sm_comm);
    v3 = dimCube(nzdim,(nydim+1),nxdim,bc_width, &sm_win_v3, &sm_comm);

    MPI_Win sm_win_w1;
    MPI_Win sm_win_w2;
    MPI_Win sm_win_w3;

    w1 = dimCube((nzdim+1),nydim,nxdim,bc_width, &sm_win_w1, &sm_comm);
    w2 = dimCube((nzdim+1),nydim,nxdim,bc_width, &sm_win_w2, &sm_comm);
    w3 = dimCube((nzdim+1),nydim,nxdim,bc_width, &sm_win_w3, &sm_comm);


    MPI_Win sm_win_tplot,sm_win_ro_u,sm_win_ro_w;
    if (myRank == root) {
        MPI_Win_allocate_shared((MPI_Aint) nz * ny * nx * sizeof(float),sizeof(float),MPI_INFO_NULL,sm_comm,&tplot,&sm_win_tplot);
        MPI_Win_allocate_shared((MPI_Aint) nzdim   *      sizeof(real), sizeof(real), MPI_INFO_NULL,sm_comm,&ro_u, &sm_win_ro_u);
        MPI_Win_allocate_shared((MPI_Aint) (nzdim+1) *    sizeof(real), sizeof(real), MPI_INFO_NULL,sm_comm,&ro_w, &sm_win_ro_w);
    } else {
        MPI_Win_allocate_shared((MPI_Aint) 0                           ,sizeof(float),MPI_INFO_NULL,sm_comm,&tplot,&sm_win_tplot);
        MPI_Win_allocate_shared((MPI_Aint) 0                           ,sizeof(real), MPI_INFO_NULL,sm_comm,&ro_u, &sm_win_ro_u);
        MPI_Win_allocate_shared((MPI_Aint) 0                           ,sizeof(real), MPI_INFO_NULL,sm_comm,&ro_w, &sm_win_ro_w);
    } // end if //

    MPI_Aint sz;
    int dispUnit;

    MPI_Win_shared_query(sm_win_tplot, MPI_PROC_NULL, &sz,&dispUnit,&tplot);
    MPI_Win_shared_query(sm_win_ro_u,  MPI_PROC_NULL, &sz,&dispUnit,&ro_u);
    MPI_Win_shared_query(sm_win_ro_w,  MPI_PROC_NULL, &sz,&dispUnit,&ro_w);
    // end of allocating memory for arrays and matrices

    MPI_Barrier(sm_comm);
    elapsed_time = -MPI_Wtime();


    /*
     * Set and plot the initial condition
     */
    MPI_Win_lock_all(0,sm_win_t1);
    MPI_Win_lock_all(0,sm_win_t2);
    MPI_Win_lock_all(0,sm_win_tplot);
    MPI_Win_lock_all(0,sm_win_p1);
    MPI_Win_lock_all(0,sm_win_p3);
    MPI_Win_lock_all(0,sm_win_u1);
    MPI_Win_lock_all(0,sm_win_u3);
    MPI_Win_lock_all(0,sm_win_v1);
    MPI_Win_lock_all(0,sm_win_v3);
    MPI_Win_lock_all(0,sm_win_w1);
    MPI_Win_lock_all(0,sm_win_w3);
    MPI_Win_lock_all(0,sm_win_ro_u);
    MPI_Win_lock_all(0,sm_win_ro_w);


    if (myRank == root) {
        ic(t1,u1,v1,w1,ro_u,ro_w,dx,dy,dz,deltau,i1,i2,j1,j2,k1,k2,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);
    } // end if //


    MPI_Win_sync(sm_win_t1);
    MPI_Win_sync(sm_win_u1);
    MPI_Win_sync(sm_win_v1);
    MPI_Win_sync(sm_win_w1);
    MPI_Win_sync(sm_win_ro_u);
    MPI_Win_sync(sm_win_ro_w);
    MPI_Barrier(sm_comm);


    //#pragma omp parallel

    //nthreads = omp_get_num_threads();
    bc(u1,v1,w1,i1,i2,j1,j2,k1,k2,bc_width,myRank,commSize);
    bc4T(t1,t2,i1,i2,j1,j2,k1,k2,bc_width,myRank,commSize);

    MPI_Win_sync(sm_win_u1);
    MPI_Win_sync(sm_win_v1);
    MPI_Win_sync(sm_win_w1);
    MPI_Win_sync(sm_win_t1);
    MPI_Win_sync(sm_win_t2);
    MPI_Barrier(sm_comm);


    const int lk1 = BLOCK_LOW (myRank,commSize,nzdim);
    const int lk2 = BLOCK_HIGH(myRank,commSize,nzdim)+1;
    //printf("myRank: %3d, lk1: %3d, lk2: %3d, ndimz: %3d\n", myRank, lk1,lk2, nzdim);

    // acctually only the BC of t1 should be copyed to t2
    //#pragma omp for
    for (int level=lk1; level<lk2; ++level) {
        for (int row=0; row<nydim; ++row) {
            for (int col=0; col<nxdim; ++col) {
                t2[level][row][col]=t1[level][row][col];
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_t2);
    MPI_Barrier(sm_comm);


    const int ii1 = BLOCK_LOW (myRank,commSize,(1+i2-i1)) + bc_width;
    const int ii2 = BLOCK_HIGH(myRank,commSize,(1+i2-i1)) + bc_width;
    //printf("myRank: %3d, --> ii1: %3d, ii2: %3d, i1: %3d, i2: %3d\n", myRank,ii1,ii2, i1,i2);


    // Ploting initial conditions //

    //#pragma omp for
    for (int col=ii1; col<=ii2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (t1[level][row][col] - thetaBar);
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_tplot);
    MPI_Barrier(sm_comm);
    //#pragma omp single
    if (myRank == root) putfield("T",0.0,tplot, nx,ny,nz);
    MPI_Barrier(sm_comm);


    //#pragma omp for
    for (int col=ii1; col<=ii2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(u1[level][row][col]+u1[level][row][col+1]));
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_tplot);
    MPI_Barrier(sm_comm);
    //#pragma omp single
    if (myRank == root) putfield("U",0.0,tplot, nx,ny,nz);
    MPI_Barrier(sm_comm);


    //#pragma omp for
    for (int col=ii1; col<=ii2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(v1[level][row][col]+v1[level][row+1][col]));
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_tplot);
    MPI_Barrier(sm_comm);
    //#pragma omp single
    if (myRank == root) putfield("V",0.0,tplot, nx,ny,nz);
    MPI_Barrier(sm_comm);

    //#pragma omp for
    for (int col=ii1; col<=ii2; ++col) {
        for (int row=j1; row<=j2; ++row) {
            for (int level=k1; level<=k2; ++level) {
                tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(w1[level][row][col]+w1[level+1][row][col]));
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_tplot);
    MPI_Barrier(sm_comm);
    //#pragma omp single
    if (myRank == root) putfield("W",0.0,tplot, nx,ny,nz);
    MPI_Barrier(sm_comm);


    const int lS = BLOCK_LOW (myRank,commSize,nx);
    const int lE = BLOCK_HIGH(myRank,commSize,nx) + 1;
    //printf("myRank: %3d, lS: %3d, lE: %3d, nx: %3d\n", myRank,lS,lE, nx);

    //#pragma omp for
    for (int col=lS; col<lE; ++col) {
        for (int row=0; row<ny; ++row) {
            for (int level=0; level<nz; ++level) {
                tplot[col*nz*ny + row*nz + level] = (float) (p1[level][row][col]);
            } // end for //
        } // end for //
    } // end for //
    MPI_Win_sync(sm_win_tplot);
    MPI_Barrier(sm_comm);
    //#pragma omp single
    if (myRank == root) putfield("P",0.0,tplot, nx,ny,nz);
    MPI_Barrier(sm_comm);

    // end of Plotting initial conditions //


    /*
     * .. Integrate .....................................................
     */

    const int lk3 = BLOCK_LOW (myRank,commSize,(nzdim+1));
    const int lk4 = BLOCK_HIGH(myRank,commSize,(nzdim+1));
    //printf("myRank: %3d, lk3: %3d, lk4: %3d, ndimz: %3d\n", myRank, lk3,lk4, nzdim);

    const int kk1 = BLOCK_LOW (myRank,commSize,(1+k2-k1)) + bc_width;
    const int kk2 = BLOCK_HIGH(myRank,commSize,(1+k2-k1)) + bc_width;
    //printf("myRank: %3d, --> kk1: %3d, kk2: %3d, k1: %3d, k2: %3d\n", myRank,kk1,kk2, k1,k2);

    tstep=dt;
    for (int n=1 ; n<=nstep; n++) {

        //#pragma omp for nowait
        for (int level=lk1; level<lk2; ++level) {
            for (int row=0; row<nydim; ++row) {
                for (int col=0; col<=nxdim; ++col) {
                    u3[level][row][col]=u1[level][row][col];
                } // end for //
            } // end for //
        } // end for //

        //#pragma omp for nowait
        for (int level=lk1; level<lk2; ++level) {
            for (int row=0; row<=nydim; ++row) {
                for (int col=0; col<nxdim; ++col) {
                    v3[level][row][col]=v1[level][row][col];
                } // end for //
            } // end for //
        } // end for //



        //#pragma omp for // using the implicit barrier to syncromize
        for (int level=lk3; level<=lk4; ++level) {
            for (int row=0; row<nydim; ++row) {
                for (int col=0; col<nxdim; ++col) {
                    w3[level][row][col]=w1[level][row][col];
                } // end for //
            } // end for //
        } // end for //

        MPI_Win_sync(sm_win_u3);
        MPI_Win_sync(sm_win_v3);
        MPI_Win_sync(sm_win_w3);
        MPI_Barrier(sm_comm);

        //  . . . Compute values at next step                           //

        advection(t2,u3,v3,w3,t1,u2,v2,w2,dx,dy,dz,dt,i1,i2,j1,j2,k1,k2,nxdim,nydim,nzdim,tstep,advection_type,sm_comm,&sm_win_t2);
        //MPI_Win_sync(sm_win_t2);  /// OOOOOJJJJOOOO///
        MPI_Win_sync(sm_win_u3);
        MPI_Win_sync(sm_win_v3);
        MPI_Win_sync(sm_win_w3);
        MPI_Barrier(sm_comm);



        diffusion(t2,u3,v3,w3,t1,u1,v1,w1,dx*dx,dy*dy,dz*dz,i1,i2,j1,j2,kk1,kk2,dt,tstep,ku,kv,kw,kt);// KU,KV,KW,KT);
        MPI_Win_sync(sm_win_t2);
        MPI_Win_sync(sm_win_u3);
        MPI_Win_sync(sm_win_v3);
        MPI_Win_sync(sm_win_w3);
        MPI_Barrier(sm_comm);

        pgf(p3,u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,i1,i2,j1,j2,k1,k2,tstep, g,thetaBar,cs*cs,bc_width,sm_comm,&sm_win_u3,&sm_win_v3,&sm_win_w3);

        //  . . . Set boundary conditions for T                     //
        bc4T(t2,t1,i1,i2,j1,j2,k1,k2,bc_width,myRank,commSize);
        MPI_Win_sync(sm_win_p3);
        MPI_Win_sync(sm_win_t1);
        MPI_Win_sync(sm_win_t2);
        MPI_Barrier(sm_comm);

        if (n%nplot == 0) {
            //#pragma omp for
            for (int col=ii1; col<=ii2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (t2[level][row][col] - thetaBar);
                    } // end for //
                } // end for //
            } // end for //
            MPI_Win_sync(sm_win_tplot);
            MPI_Barrier(sm_comm);
            //#pragma omp single
            if (myRank == root) putfield("T",(float) dt*n,tplot, nx,ny,nz);
            MPI_Barrier(sm_comm);



            //#pragma omp for
            for (int col=ii1; col<=ii2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(u3[level][row][col]+u3[level][row][col+1]));
                    } // end for //
                } // end for //
            } // end for //
            MPI_Win_sync(sm_win_tplot);
            MPI_Barrier(sm_comm);
            //#pragma omp single
            if (myRank == root) putfield("U",(float) dt*n,tplot, nx,ny,nz);
            MPI_Barrier(sm_comm);

            //#pragma omp for
            for (int col=ii1; col<=ii2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(v3[level][row][col]+v3[level][row+1][col]));
                    } // end for //
                } // end for //
            } // end for //
            MPI_Win_sync(sm_win_tplot);
            MPI_Barrier(sm_comm);
            //#pragma omp single
            if (myRank == root) putfield("V",(float) dt*n,tplot, nx,ny,nz);
            MPI_Barrier(sm_comm);

            //#pragma omp for
            for (int col=ii1; col<=ii2; ++col) {
                for (int row=j1; row<=j2; ++row) {
                    for (int level=k1; level<=k2; ++level) {
                        tplot[(col-i1)*nz*ny + (row-j1)*nz + (level-k1)] = (float) (0.5*(w3[level][row][col]+w3[level+1][row][col]));
                    } // end for //
                } // end for //
            } // end for //
            MPI_Win_sync(sm_win_tplot);
            MPI_Barrier(sm_comm);
            //#pragma omp single
            if (myRank == root) putfield("W",(float) dt*n,tplot, nx,ny,nz);
            MPI_Barrier(sm_comm);

            //#pragma omp for
            for (int col=lS; col<lE; ++col) {
                for (int row=0; row<ny; ++row) {
                    for (int level=0; level<nz; ++level) {
                        tplot[col*nz*ny + row*nz + level] = (float) (p3[level][row][col]);
                    } // end for //
                } // end for //
            } // end for //
            MPI_Win_sync(sm_win_tplot);
            MPI_Barrier(sm_comm);
            //#pragma omp single
            if (myRank == root) putfield("P",(float) dt*n,tplot, nx,ny,nz);
            MPI_Barrier(sm_comm);

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


        MPI_Win_sync(sm_win_t1);
        MPI_Win_sync(sm_win_t2);
        MPI_Win_sync(sm_win_u3);
        MPI_Win_sync(sm_win_u1);
        MPI_Win_sync(sm_win_v3);
        MPI_Win_sync(sm_win_v1);
        MPI_Win_sync(sm_win_w3);
        MPI_Win_sync(sm_win_w1);
        MPI_Win_sync(sm_win_p3);
        MPI_Win_sync(sm_win_p1);
        MPI_Barrier(sm_comm);
    }	// end of time loop n = 1,...,nstep //


    MPI_Win_unlock_all(sm_win_t1);
    MPI_Win_unlock_all(sm_win_t2);
    MPI_Win_unlock_all(sm_win_tplot);
    MPI_Win_unlock_all(sm_win_p1);
    MPI_Win_unlock_all(sm_win_p3);
    MPI_Win_unlock_all(sm_win_u1);
    MPI_Win_unlock_all(sm_win_u3);
    MPI_Win_unlock_all(sm_win_v1);
    MPI_Win_unlock_all(sm_win_v3);
    MPI_Win_unlock_all(sm_win_w1);
    MPI_Win_unlock_all(sm_win_w3);
    MPI_Win_unlock_all(sm_win_ro_u);
    MPI_Win_unlock_all(sm_win_ro_w);

    MPI_Barrier(sm_comm);
    elapsed_time += MPI_Wtime();

    if (myRank == root) {
        printf ("\n\nIt tooks %14.6e seconds for %d processes to finish\n", elapsed_time, commSize);
        if (sizeof(real) == 8) {
            printf("Double precision version\n");
        } else {
            printf("Single precision version\n");
        } // end if
    } // end if

    freeCube(&t2, &sm_win_t2);
    freeCube(&t1, &sm_win_t1);

    freeCube(&p1, &sm_win_p1);
    freeCube(&p2, &sm_win_p2);
    freeCube(&p3, &sm_win_p3);

    freeCube(&u1, &sm_win_u1);
    freeCube(&u2, &sm_win_u2);
    freeCube(&u3, &sm_win_u3);

	freeCube(&v1, &sm_win_v1);
	freeCube(&v2, &sm_win_v2);
	freeCube(&v3, &sm_win_v3);

    freeCube(&w1, &sm_win_w1);
    freeCube(&w2, &sm_win_w2);
    freeCube(&w3, &sm_win_w3);

    MPI_Win_free(&sm_win_tplot);
    MPI_Win_free(&sm_win_ro_u);
    MPI_Win_free(&sm_win_ro_w);
    MPI_Finalize();


    return 0;
} // end main() //

