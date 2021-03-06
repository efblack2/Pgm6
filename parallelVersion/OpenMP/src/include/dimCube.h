#include <omp.h>
#include "real.h"

real*** dimCube(int high, int row, int col)
{

    real ***cube;

    cube       = (real ***)  malloc(  high  * sizeof(real**));
    cube[0]    = (real  **)  malloc(  high  * row  * sizeof(real*));
    cube[0][0] = (real   *)  malloc(( high  * row * col) * sizeof(real));

    // creating the cube array
    for(int k=0; k < high; ++k){
        cube[k] = cube[0]  + (k * row);
        for(int j=0; j < row; ++j){
            cube[k][j] = cube[0][0]  +  (k *  row * col )  + j * col;
        } // end for //
    } // end for //
    // end of creating the cube array

    #pragma omp parallel for
    for(int l = 0; l < high; ++l) {
        for(int r = 0; r < row ; ++r) {
            for(int c = 0; c < col; ++c) {
                cube[l][r][c] = (real) 0.0;
            } // end for //
        } // end for //
    } // end for //


   return cube;

} // end of dimCube() //

void freeCube(real ****cube)
{
    free((*cube)[0][0]);
    free((*cube)[0]);
    free((*cube));
} // end of freeCube //
