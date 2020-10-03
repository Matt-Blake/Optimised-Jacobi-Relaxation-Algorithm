/* ****************************************************************
 * poisson_test.cpp
 *
 * Part of ENCE464 Assignment 2
 * This program attemps to solve Poisson's equation for a 3-D
 * cube using Jacobi relaxation.
 *
 * Based off poisson_test.cpp - Michael Hayes, UC ECE
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 03/10/2020
 *
 * ***************************************************************/

#include <cstdio>
#include <cstdlib>
#include "poisson.hpp"

#define VOXEL_SPACING       0.1     // The spacing (in all directions) between voxels (meters)



int main (int argc, char *argv[])
{
    double *source;
    double *potential;
    unsigned int N;
    unsigned int num_iters;
    unsigned int num_cores;    
    unsigned int x_size;
    unsigned int y_size;
    unsigned int z_size;    
    double delta = VOXEL_SPACING;

    if (argc < 3)
    {
        fprintf (stderr, "Usage: %s size numiters\n", argv[0]);
        return 1;
    }

    // Set the size of the 3-D array (representing volume)  
    N = atoi(argv[1]);
    x_size = N;
    y_size = N;
    z_size = N;

    num_iters = atoi(argv[2]); // The number of iterations of Poisson's equation to perform

    if (argc > 3)
        num_cores = atoi(argv[3]); // If the user chooses to define the number of cores used
    else
        num_cores = 0;

    source = (double*) calloc(x_size * y_size * z_size, sizeof(*source));
    potential = (double*) calloc(x_size * y_size * z_size, sizeof(*potential));

    source[((z_size / 2 * y_size) + y_size / 2) * x_size + x_size / 2] = 1.0;    
    
    poisson_dirichlet(source, potential, 0, x_size, y_size, z_size, delta,
                      num_iters, num_cores);

    return 0;
}
