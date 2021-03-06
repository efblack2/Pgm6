
#include <mpi.h>
#include "myMPI.h"
#include "real.h"

real*** dimCube(int level, int row, int col, int bw, MPI_Win *sm_win, MPI_Comm *sm_comm)
{

    int myRank,commSize;
    MPI_Comm_rank(*sm_comm,&myRank);
    MPI_Comm_size(*sm_comm,&commSize);

    const int size = level  * row * col;
    real ***cube;

    int zdim = level-2*bw;

    int low = BLOCK_LOW (myRank,commSize,zdim)+bw;
    int hi  = BLOCK_HIGH(myRank,commSize,zdim)+bw;

    if (myRank == 0) low-=bw;
    if (myRank == commSize-1) hi+=bw;



    cube       = (real ***)  malloc(  level         * sizeof(real**));
    cube[0]    = (real  **)  malloc(  level  * row  * sizeof(real*));
    //cube[0][0] = (real   *)  calloc(( high  * row * col),sizeof(real));

    if (myRank == 0) {
        MPI_Win_allocate_shared((MPI_Aint) size*sizeof(real), sizeof(real), MPI_INFO_NULL,*sm_comm,&cube[0][0],sm_win);
    } else {
        MPI_Win_allocate_shared((MPI_Aint) 0,                 sizeof(real), MPI_INFO_NULL,*sm_comm,&cube[0][0],sm_win);
    } // end if //



    MPI_Barrier(*sm_comm);  // is this really needed? //

    MPI_Aint sz;
    int dispUnit;
    MPI_Win_shared_query(*sm_win, MPI_PROC_NULL, &sz,&dispUnit,&cube[0][0]);

    for(int l=0; l < level; ++l){
        cube[l] = cube[0]  + (l * row);
        for(int r=0; r < row; ++r){
            cube[l][r] = cube[0][0]  +  (l *  row * col )  + r * col;
        } // end for //
    } // end for //

    MPI_Barrier(*sm_comm);  // is this really needed? //


    MPI_Win_lock_all(0,*sm_win);

    for (int l=low; l<hi; ++l){
        for (int r=0; r<row; ++r){
            for (int c=0; c<col; ++c){
                cube[l][r][c] = (real) 0.0;
            } // end for //
        } // end for //
    } // end for //

    MPI_Win_sync(*sm_win);
    MPI_Barrier(*sm_comm);
    MPI_Win_unlock_all(*sm_win);

    MPI_Barrier(*sm_comm);  // is this really needed? //

   return cube;

} // end of dimCube() //

void freeCube(real ****cube, MPI_Win *sm_win)
{
    MPI_Win_free(sm_win);
    //free((*cube)[0][0]);
    free((*cube)[0]);
    free((*cube));
} // end of freeCube //
