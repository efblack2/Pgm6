#include <stdlib.h>
#include "real.h"
/*
 * ========================= update ====================
 * Update: replace old values with new ones
 * We are NOT copying ghost points here.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1,q2	real arrays	old, new data arrays
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		size of arrays
 *
 */

void update(real ****cube1, real ****cube2)
{

    real ***temp;
    
    temp = *cube1;
    *cube1 = *cube2;
    *cube2 = temp;
    
} // end update() //

