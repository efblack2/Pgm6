/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm1:  Linear and nonlinear advection
 *  >>>>> PUT YOUR NAME HERE!  Yes, you!! <<<<<
 */

#define MAIN
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <sys/time.h>
#include <cuda.h> 

#include "std2file.hpp"
#include "prototypes.cuh"
#include "real.h"
#include "dimCube.h"
#include "fatal.hpp"
#include "constant.cuh"

/*

    timeGPU=0.0;
    istat = cudaEventRecord(startEvent,0);   

    istat = cudaEventRecord(stopEvent,0);
    istat = cudaEventSynchronize(stopEvent);
    istat = cudaEventElapsedTime(&tempGPU, startEvent, stopEvent);
    timeGPU+=tempGPU;

*/

/*
 * Definitions
 * NX is the number of physical grid points - Not including 'ghost points'
 * bc_width is the number of ghost points - on each side of the physical grid
 * iW is the first C index to use - the left-most physical point
 * iE is the last C index to use - the right-most physical point
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
    // Input and output streams
    ifstream infile;
    ofstream outfile;
    // end of Input and output streams

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
    real ***p;
    real ***u;
    real ***v;
    real ***w;
	real cs,deltau;
	real KU,KV,KW,KTemp;
    int nx,ny,nz,iW, iE, jS,jN, kB, kT, nxdim,nydim, nzdim, bc_width=0;
    int nstep,nplot;
    char advection_type;

// Cuda declarations
    //cudaEvent_t startEvent,stopEvent;
    //float tempGPU, timeGPU;
    //cudaDeviceProp prop;
    cudaError_t istat;
    dim3 grid,block;
    
    // Defining block sizes
    block.x = TILESIZEX;
    block.y = TILESIZEY;
    block.z = 1;
    // end of defining block sizes
    
    //dim3 grid[3], block[3];
    real *ro_u_d, *ro_w_d;
    real *t1_d, *t2_d, *tTemp_d;
    real *p1_d, *p2_d, *p3_d;
    real *u1_d, *u2_d, *u3_d;
    real *v1_d, *v2_d, *v3_d;
    real *w1_d, *w2_d, *w3_d;
    //istat = cudaEventCreate(&startEvent);
    //istat = cudaEventCreate(&stopEvent);
// end of Cuda declarations



    /* Parameters and input .................................... */

    if (argc > 1 ) {
        infile.open(argv[1], ios::in);
        if (!infile) {
            cerr <<" Can't open input file " << argv[1] << '!' << endl;
            exit(-1);
        } // end if // 
        std2file(cin,infile);
    } else {
        cout << "Use: " << argv[0] << " filename \n";
        cout << "Example: " << argv[0] << "  data.txt \n"; 
        cout << "advection_scheme can be: L for Lax-Wendroff, C for Crowley, T for Takacs, P for Piecewise, U for Upstream\n";
        exit(-1);
    } // endif //
    
    if (sizeof(real) == 8) {
        cout << "Double precision version\n";
    } else {
        cout << "Single precision version\n";
    } // end if


    cin >> nstep >> nplot >> advection_type;
    cin >> nx >> ny >> nz;
    cin >> dx >> dy >> dz >> dt;
    cin >> cs;
    cin >> KU >> KV >> KW >> KTemp;
    cin >> deltau;
    cin >> x0[0] >> y0[0] >> z0[0];
    cin >> rx[0] >> ry[0] >> rz[0];
    cin >> deltaTheta[0] >>  deltaV[0];
    cin >> x0[1] >> y0[1] >> z0[1];
    cin >> rx[1] >> ry[1] >> rz[1];
    cin >> deltaTheta[1] >>  deltaV[1];
    


    cout << "Program 6, ATMS 502/CSE 566, Fall 2016\n";
    cout << "Domain (" << setw(3) << nx << ","
                       << setw(3) << ny << ","     
                       << setw(3) << nz << ")"  << endl;

    cout << "dx=" << fixed << setprecision(5) << dx << endl;
    cout << "dy=" << fixed << setprecision(5) << dy << endl;
    cout << "dz=" << fixed << setprecision(5) << dz << endl; 
    cout << "dt=" << fixed << setprecision(5) << dt << endl;
    cout << "cs=" << fixed << setprecision(5) << cs << endl;  
    cout << "KU=" << fixed << setprecision(5) << KU << endl;
    cout << "KV=" << fixed << setprecision(5) << KV << endl;
    cout << "KW=" << fixed << setprecision(5) << KW << endl;
    cout << "KT=" << fixed << setprecision(5) << KTemp << endl;
    cout << "deltau=" << fixed << setprecision(5) << deltau << endl;
                       
	
    cout << "x_0=" << fixed << setprecision(5) << x0[0] << endl;
    cout << "y_0=" << fixed << setprecision(5) << y0[0] << endl;
    cout << "z_0=" << fixed << setprecision(5) << z0[0] << endl;
    cout << "rx_0=" << fixed << setprecision(5) << rx[0] << endl;
    cout << "ry_0=" << fixed << setprecision(5) << ry[0] << endl;
    cout << "rz_0=" << fixed << setprecision(5) << rz[0] << endl;
    cout << "deltaTheta_0=" << fixed << setprecision(5) << deltaTheta[0] << endl;
    cout << "deltaV_0=" << fixed << setprecision(5) << deltaV[0] << endl;

    cout << "x_1=" << fixed << setprecision(5) << x0[1] << endl;
    cout << "y_1=" << fixed << setprecision(5) << y0[1] << endl;
    cout << "z_1=" << fixed << setprecision(5) << z0[1] << endl;
    cout << "rx_1=" << fixed << setprecision(5) << rx[1] << endl;
    cout << "ry_1=" << fixed << setprecision(5) << ry[1] << endl;
    cout << "rz_1=" << fixed << setprecision(5) << rz[1] << endl;
    cout << "deltaTheta_1=" << fixed << setprecision(5) << deltaTheta[1] << endl;
    cout << "deltaV_1=" << fixed << setprecision(5) << deltaV[1] << endl;

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
        cout << "Integrate: Error, unrecognized advection type " << advection_type << endl;
        exit(-1);
        break;
    } // end switch //
    iW          = bc_width;
    iE          = iW+nx-1;
    jS          = bc_width;
    jN          = jS+ny-1;
    kB          = bc_width;
    kT          = kB+nz-1;
    nxdim       = nx+2*bc_width;
    nydim       = ny+2*bc_width;
    nzdim       = nz+2*bc_width;

    // copying constant memory to device
    istat = cudaMemcpyToSymbol(::nx,&nx, sizeof(int));
    istat = cudaMemcpyToSymbol(::ny,&ny, sizeof(int));
    istat = cudaMemcpyToSymbol(::nz,&nz, sizeof(int));
    istat = cudaMemcpyToSymbol(::bcw,&bc_width, sizeof(int));
    istat = cudaMemcpyToSymbol(::dx,&dx, sizeof(real));
    istat = cudaMemcpyToSymbol(::dy,&dy, sizeof(real));
    istat = cudaMemcpyToSymbol(::dz,&dz, sizeof(real));
    istat = cudaMemcpyToSymbol(::dt,&dt, sizeof(real));
    cs=cs*cs;
    istat = cudaMemcpyToSymbol(cs2,&cs, sizeof(real));
    istat = cudaMemcpyToSymbol(::KU,&KU, sizeof(real));
    istat = cudaMemcpyToSymbol(::KV,&KV, sizeof(real));
    istat = cudaMemcpyToSymbol(::KW,&KW, sizeof(real));
    istat = cudaMemcpyToSymbol(::KTemp,&KTemp, sizeof(real));
    istat = cudaMemcpyToSymbol(::iW,&iW, sizeof(int));
    istat = cudaMemcpyToSymbol(::iE,&iE, sizeof(int));
    istat = cudaMemcpyToSymbol(::jS,&jS, sizeof(int));
    istat = cudaMemcpyToSymbol(::jN,&jN, sizeof(int));
    istat = cudaMemcpyToSymbol(::kB,&kB, sizeof(int));
    istat = cudaMemcpyToSymbol(::kT,&kT, sizeof(int));
    istat = cudaMemcpyToSymbol(::nxdim,&nxdim, sizeof(int));
    istat = cudaMemcpyToSymbol(::nydim,&nydim, sizeof(int));
    istat = cudaMemcpyToSymbol(::nzdim,&nzdim, sizeof(int));
    istat = cudaMemcpyToSymbol(::g,&g, sizeof(real));
    istat = cudaMemcpyToSymbol(::thetaBar,&thetaBar, sizeof(real));
    istat = cudaMemcpyToSymbol(::advection_type,&advection_type, sizeof(char));
    
    if(istat != cudaSuccess) FATAL("Unable to copy constant memory to device");
    // end of copying constant memory to device

    
    
    t1 = dimCube(nzdim,nydim,nxdim);
    t2 = dimCube(nzdim,nydim,nxdim);
    istat = cudaMalloc((void **) (size_t) &t1_d, (size_t) nzdim*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for t1_d into device");
    istat = cudaMalloc((void **) (size_t) &t2_d, (size_t) nzdim*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for t2_d into device");
    istat = cudaMalloc((void **) (size_t) &tTemp_d, (size_t) nzdim*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for tTemp_d into device");
    
	tplot = (float *) malloc( nx * ny * nz * sizeof(float*));

    p = dimCube(nz,ny,nx);

    istat = cudaMalloc((void **) &p1_d, (size_t) nz*ny*nx*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for p1_d into device");

    istat = cudaMalloc((void **) &p2_d, (size_t) nz*ny*nx*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for p2_d into device");
    
    istat = cudaMalloc((void **) &p3_d, (size_t) nz*ny*nx*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for p3_d into device");

    u = dimCube(nzdim,nydim,nxdim+1);

    istat = cudaMalloc((void **) &u1_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for u1_d into device");

    istat = cudaMalloc((void **) &u2_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for u2_d into device");

    istat = cudaMalloc((void **) &u3_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for u3_d into device");
    
    v = dimCube(nzdim,nydim+1,nxdim);
    
    istat = cudaMalloc((void **) &v1_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for v1_d into device");

    istat = cudaMalloc((void **) &v2_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for v2_d into device");

    istat = cudaMalloc((void **) &v3_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for v3_d into device");

    w = dimCube(nzdim+1,nydim,nxdim);
    
    istat = cudaMalloc((void **) &w1_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for w1_d into device");

    istat = cudaMalloc((void **) &w2_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for w2_d into device");

    istat = cudaMalloc((void **) &w3_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for w3_d into device");


    ro_u = (real *) malloc(nzdim*sizeof(real));
    ro_w = (real *) malloc((nzdim+1)*sizeof(real));

    istat = cudaMalloc((void **) &ro_u_d, (size_t) nzdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for ro_u_d into device");

    istat = cudaMalloc((void **) &ro_w_d, (size_t) (nzdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to allocate memory for ro_w_d into device");
    
    
    
    // end of allocating memory for arrays and matrices


    gettimeofday(&tp,NULL);
    elapsed_time = -(tp.tv_sec*1.0e6 + tp.tv_usec);  

    /*
     * Set and plot the initial condition
     */
     
    ic(t1,u,v,w,ro_u,ro_w,dx,dy,dz,deltau,iW,iE,jS,jN,kB,kT,bc_width,nx,ny,nz,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);


    istat = cudaMemcpy(t1_d, &t1[0][0][0],(size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for t1_d into device");

    istat = cudaMemcpy(t2_d, &t2[0][0][0],(size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for t2_d into device");

    istat = cudaMemset(tTemp_d, 0 ,(size_t) nzdim*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to set memory for tTemp_d into device");


    istat = cudaMemcpy(p1_d, &p[0][0][0],(size_t) nz*ny*nx*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for p1_d into device");


    istat = cudaMemset(p2_d, 0,(size_t) nz*ny*nx*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for p2_d into device");

    istat = cudaMemset(p3_d, 0,(size_t) nz*ny*nx*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for p3_d into device");


    istat = cudaMemcpy(u1_d, &u[0][0][0],(size_t) nzdim*nydim*(nxdim+1)*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for u1_d into device");
    
    istat = cudaMemset(u2_d, 0,(size_t) nzdim*nydim*(nxdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for u2_d into device");
    
    istat = cudaMemset(u3_d, 0,(size_t) nzdim*nydim*(nxdim+1)*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for u3_d into device");
    


    istat = cudaMemcpy(v1_d, &v[0][0][0],(size_t) nzdim*(nydim+1)*nxdim*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for v1_d into device");

    istat = cudaMemset(v2_d, 0,(size_t) nzdim*(nydim+1)*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for v2_d into device");

    istat = cudaMemset(v3_d, 0,(size_t) nzdim*(nydim+1)*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for v3_d into device");



    istat = cudaMemcpy(w1_d, &w[0][0][0],(size_t) (nzdim+1)*nydim*nxdim*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for w1_d into device");

    istat = cudaMemset(w2_d, 0,(size_t) (nzdim+1)*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for w2_d into device");

    istat = cudaMemset(w3_d, 0,(size_t) (nzdim+1)*nydim*nxdim*sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to clear memory for w3_d into device");


    istat = cudaMemcpy(ro_u_d, ro_u,(size_t) nzdim*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for ro_u_d into device");

    istat = cudaMemcpy(ro_w_d, ro_w,(size_t) (nzdim+1)*sizeof(real), cudaMemcpyHostToDevice);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for ro_w_d into device");



    // Defining grid and block sizes
    grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
    grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
    // end of defining grid and block sizes        
    bcTB<<<grid,block>>>(u1_d,v1_d,w1_d,t1_d,t2_d);


    // Defining grid and block sizes
    grid.x  = (ny + TILESIZEX-1) / TILESIZEX;
    grid.y  = (nz + TILESIZEY-1) / TILESIZEY;
    // end of defining grid and block sizes        
    bcEW<<<grid,block>>>(u1_d,v1_d,w1_d,t1_d,t2_d);


    // Defining grid and block sizes
    grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
    //grid.y  = (nz + TILESIZEY-1) / TILESIZEY;
    // end of defining grid and block sizes        
    bcNS<<<grid,block>>>(u1_d,v1_d,w1_d,t1_d,t2_d);





    /* eliminar luego */
    istat = cudaMemcpy(&u[0][0][0],u1_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real), cudaMemcpyDeviceToHost);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for u1_d back to host");

    istat = cudaMemcpy(&v[0][0][0],v1_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for v1_d back to host");

    istat = cudaMemcpy(&w[0][0][0],w1_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for w1_d back to host");

    istat = cudaMemcpy(&t1[0][0][0],t1_d, (size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for t1_d back to host");

    istat = cudaMemcpy(&t2[0][0][0],t2_d, (size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
    if(istat != cudaSuccess) FATAL("Unable to copy memory for t2_d back to host");
    /* eliminar luego */



    
    
    // acctually only the BC of t1 should be copyed to t2 
    
    for (int level=0; level<nzdim; ++level) {
        for (int row=0; row<nydim; ++row) {
            for (int col=0; col<nxdim; ++col) {
                t2[level][row][col]=t1[level][row][col];
            } // end for //
        } // end for //
    } // end for //
    
    // Ploting initial conditions //
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (t1[level][row][col] - thetaBar);
            } // end for //
        } // end for //
    } // end for // 
   
    putfield((char *) "T",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(u[level][row][col]+u[level][row][col+1]));
            } // end for //
        } // end for //
    } // end for // 
   
    putfield((char *) "U",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(v[level][row][col]+v[level][row+1][col]));
            } // end for //
        } // end for //
    } // end for //     
    
    putfield((char *) "V",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(w[level][row][col]+w[level+1][row][col]));
            } // end for //
        } // end for //
    } // end for //
     
    putfield((char *) "W",0.0,tplot, nx,ny,nz);

   
    for (int col=0; col<nx; ++col) {
        for (int row=0; row<ny; ++row) {
            for (int level=0; level<nz; ++level) {
                tplot[col*nz*ny + row*nz + level] = (float) (p[level][row][col]);
            } // end for //
        } // end for //
    } // end for //
    
    putfield((char *) "P",0.0,tplot, nx,ny,nz);

    // end of Plotting initial conditions //


    /*
     * .. Integrate .....................................................
     */
    
    //timeGPU=0.0;
     
    tstep=dt;
    istat = cudaMemcpyToSymbol(::tstep,&tstep, sizeof(real));
    if(istat != cudaSuccess) FATAL("Unable to copy constant tstep to device");
    
    for (int n=1 ; n<=nstep; n++) {
        

        //  . . . Compute values at next step                           //
        istat = cudaMemcpy(u3_d,u1_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real), cudaMemcpyDeviceToDevice);
        if(istat != cudaSuccess) FATAL("Unable to copy update u3 in device");

        istat = cudaMemcpy(v3_d,v1_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real), cudaMemcpyDeviceToDevice);
        if(istat != cudaSuccess) FATAL("Unable to copy update v3 in device");

        istat = cudaMemcpy(w3_d,w1_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToDevice);
        if(istat != cudaSuccess) FATAL("Unable to copy update w3 in device");

        // needed here for solution proposed for advection routines
        istat = cudaMemcpy(tTemp_d,t1_d, (size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToDevice);
        if(istat != cudaSuccess) FATAL("Unable to copy memory for tTemp_d into device");


        // Defining grid and block sizes
        grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        advectionEW<<<grid,block>>>(t2_d,t1_d,u3_d);


        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        //grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        advectionNS<<<grid,block>>>(tTemp_d,t2_d,v3_d);

        
        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        grid.y  = (nz + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        advectionTB<<<grid,block>>>(t2_d,tTemp_d,w3_d);        


        
        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        advecUVW<<<grid,block>>>(u3_d,v3_d,w3_d,u2_d,v2_d,w2_d);        
        
        


        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        //grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes            
        diffusion<<<grid,block>>>(t2_d,u3_d,v3_d,w3_d,t1_d,u1_d,v1_d,w1_d);


        
        
        //pgf1(u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,iW,iE,jS,jN,kB,kT,tstep, g,thetaBar,bc_width);
        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        //grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes            
        pgf1<<<grid,block>>>(u3_d,v3_d,w3_d,t2_d,p1_d,ro_u_d,ro_w_d);


     
        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        //grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes            
        bcTB<<<grid,block>>>(u3_d,v3_d,w3_d,t2_d,t1_d);



        // Defining grid and block sizes
        grid.x  = (ny + TILESIZEX-1) / TILESIZEX;
        grid.y  = (nz + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        bcEW<<<grid,block>>>(u3_d,v3_d,w3_d,t2_d,t1_d);

        // Defining grid and block sizes
        grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        //grid.y  = (nz + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes        
        bcNS<<<grid,block>>>(u3_d,v3_d,w3_d,t2_d,t1_d);




        // Defining grid and block sizes
        //grid.x  = (nx + TILESIZEX-1) / TILESIZEX;
        grid.y  = (ny + TILESIZEY-1) / TILESIZEY;
        // end of defining grid and block sizes            
        pgf2<<<grid,block>>>(p3_d,u3_d,v3_d,w3_d,p1_d,ro_u_d,ro_w_d);





        if (n%nplot == 0) {
        
            istat = cudaMemcpy(&u[0][0][0],u3_d, (size_t) nzdim*nydim*(nxdim+1)*sizeof(real), cudaMemcpyDeviceToHost);
            if(istat != cudaSuccess) FATAL("Unable to copy memory for u3_d back to host");

            istat = cudaMemcpy(&v[0][0][0],v3_d, (size_t) nzdim*(nydim+1)*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
            if(istat != cudaSuccess) FATAL("Unable to copy memory for v3_d back to host");

            istat = cudaMemcpy(&w[0][0][0],w3_d, (size_t) (nzdim+1)*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
            if(istat != cudaSuccess) FATAL("Unable to copy memory for w3_d back to host");

            istat = cudaMemcpy(&t2[0][0][0],t2_d, (size_t) nzdim*nydim*nxdim*sizeof(real), cudaMemcpyDeviceToHost);
            if(istat != cudaSuccess) FATAL("Unable to copy memory for t2_d back to host");

            istat = cudaMemcpy(&p[0][0][0],p3_d, (size_t) nz*ny*nx*sizeof(real), cudaMemcpyDeviceToHost);
            if(istat != cudaSuccess) FATAL("Unable to copy memory for p3_d back to host");
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (t2[level][row][col] - thetaBar);
                    } // end for //
                } // end for //
            } // end for // 
            
            putfield((char *) "T",(float) dt*n,tplot, nx,ny,nz);
            
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(u[level][row][col]+u[level][row][col+1]));
                    } // end for //
                } // end for //
            } // end for // 
             
            putfield((char *) "U",(float) dt*n,tplot, nx,ny,nz);

            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(v[level][row][col]+v[level][row+1][col]));
                    } // end for //
                } // end for //
            } // end for //     
            
            putfield((char *) "V",(float) dt*n,tplot, nx,ny,nz);
            
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(w[level][row][col]+w[level+1][row][col]));
                    } // end for //
                } // end for //
            } // end for //
            
            putfield((char *) "W",(float) dt*n,tplot, nx,ny,nz);

            
            for (int col=0; col<nx; ++col) {
                for (int row=0; row<ny; ++row) {
                    for (int level=0; level<nz; ++level) {
                        tplot[col*nz*ny + row*nz + level] = (float) (p[level][row][col]);
                    } // end for //
                } // end for //
            } // end for //
            
            putfield((char *) "P",(float) dt*n,tplot, nx,ny,nz);
        
        } // end if //

        //  . . . Do array update at end of time step                   //
        if (n == 1) { 
            tstep=2.0*dt;
            istat = cudaMemcpyToSymbol(::tstep,&tstep, sizeof(real));
            if(istat != cudaSuccess) FATAL("Unable to copy constant tstep to device");            
            update(&p2_d,&p3_d);
            update(&u2_d,&u3_d);
            update(&v2_d,&v3_d);
            update(&w2_d,&w3_d);
            
        } else {
            update(&p1_d,&p2_d);
            update(&p2_d,&p3_d);
            
            update(&u1_d,&u2_d);
            update(&u2_d,&u3_d);
            
            update(&v1_d,&v2_d);
            update(&v2_d,&v3_d);
            
            update(&w1_d,&w2_d);
            update(&w2_d,&w3_d);
        } // end if //
        update(&t1_d,&t2_d);

    }	// end of time loop n = 1,...,nstep //

    //cout << "\n\ntime in advectionEW" << setw(14) << fixed << setprecision(6) << timeGPU*1.0e-3 << " seconds to finish\n";
         
        
    gettimeofday(&tp,NULL);
    elapsed_time += (tp.tv_sec*1.0e6 + tp.tv_usec);
    cout << "\n\nIt tooks " << setw(14) << scientific << setprecision(6) 
         << elapsed_time*1.0e-6 << " seconds to finish\n";

    if (sizeof(real) == 8) {
        cout << "Double precision version\n";
    } else {
        cout << "Single precision version\n";
    } // end if



    //istat = cudaEventDestroy(startEvent);
    //istat = cudaEventDestroy(stopEvent);        

    cudaFree(t1_d);
    cudaFree(t2_d);
    cudaFree(tTemp_d);

    cudaFree(p1_d);
    cudaFree(p2_d);
    cudaFree(p3_d);

    cudaFree(u1_d);
    cudaFree(u2_d);
    cudaFree(u3_d);

    cudaFree(v1_d);
    cudaFree(v2_d);
    cudaFree(v3_d);

    cudaFree(w1_d);
    cudaFree(w2_d);
    cudaFree(w3_d);

    cudaFree(ro_u_d);
    cudaFree(ro_w_d);

    
    freeCube(&t2);
    freeCube(&t1);

    freeCube(&p);
    
    freeCube(&u);

	freeCube(&v);

    freeCube(&w);

	free(tplot);
	free(ro_u);
	free(ro_w);

    
    return 0;
} // end main() //
#undef MAIN
