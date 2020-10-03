/* ****************************************************************
 * poisson.cpp
 *
 * Part of ENCE464 Assignment 2
 * Source file containing the classes and functions needed to
 * calculate the electric potential in a 3-D cube, for a given
 * charge distribution. This is solved through Poisson's equation
 * using Jacobi relaxation. 
 *
 * Based off poisson.cpp - Michael Hayes, UC ECE
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 03/10/2020
 * ***************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "poisson.hpp"

#define VOXEL_SPACING       0.1     // The spacing (in all directions) between voxels (meters)
#define V_BOUND             0       // Voltage potential on the box boundary
#define CHARGE_VALUE        1       // The value of a point charge in the charge distribution
#define ERROR_SYMBOL        1       // The value to return if an error occurs


/*
 * Function:    poisson_args_t::allocateUserInputs
 * ------------------------------
 * Allocates the user inputs to the poisson_test program. These
 * inputs are allocated to the poisson_args_t object.
 *
 * @params:
 *      - int argc: The number of arguments passed into main
 *      - char** argv: The vector of arguments passed into main.
 * 					   This should contain the 3-D size of the box,
 *                     followed by the number of Jacobi relaxation
 *                     iterations to perform, then an optional value
 *                     specifying the number of CPU cores to use.
 * @returns:
 *      - uint8_t error: A variable that will be 1 if an error has
 *                       occurred and 0 otherwise.
 */
int poisson_args_t::allocateUserInputs(int argc, char** argv)
{  
    unsigned int N;
	uint8_t error = EXIT_SUCCESS;

    // Check if enough user inputs were given
    if (argc < 3)
    {
        fprintf (stderr, "Usage: %s size numiters\n", argv[0]); // Print error
        error = ERROR_SYMBOL;

        return error;
    }

    // Set the size of the 3-D array (representing volume) 
    N = atoi(argv[1]);
    x_size = N;
    y_size = N;
    z_size = N;

    // Set the number of iterations of Jacobi relaxation to perform
    num_iters = atoi(argv[2]);

    // Allocate the number of cores based on if the user has made a definition
    if (argc > 3)
        num_cores = atoi(argv[3]); 
    else
        num_cores = 0;

    // Set values based on #defined constants
    delta = VOXEL_SPACING;
    V_bound = V_BOUND;

    return 0;
}


/*
 * Function:    poisson_args_t::allocateVolume
 * ------------------------------
 * Initalise the array's mapping the charge distribution
 * and resulting potential.
 */
void * poisson_args_t::allocateVolume(void)
{
    // Allocate 3-D arrays for the source distribution and potential
    source = (double*) calloc(x_size * y_size * z_size, sizeof(*source));
    potential = (double*) calloc(x_size * y_size * z_size, sizeof(*potential));

    // Allocate charge to the centre of the charge distribution 
    source[((z_size / 2 * y_size) + y_size / 2) * x_size + x_size / 2] = CHARGE_VALUE; 

    return NULL;
}


/*
 * Function:    poisson_args_t::initPoissonArgs
 * ------------------------------
 * Initialises a poisson_args_t object based on user inputs.
 * This object can then be solved using the poisson_dirichlet
 * member function.
 *
 * @params:
 *      - int argc: The number of arguments passed into main
 *      - char** argv: The vector of arguments passed into main.
 * 					   This should contain the 3-D size of the box,
 *                     followed by the number of Jacobi relaxation
 *                     iterations to perform, then an optional value
 *                     specifying the number of CPU cores to use.
 */
void* poisson_args_t::initPoissonArgs(int argc, char** argv)
{  
    allocateUserInputs(argc, argv); // Initalise members based on user inputs
    allocateVolume(); // Dynamically create space for the 3-D volume arrays

    return NULL;
}

/*
 * Function:    poisson_args_t::poissonDirichlet
 * ------------------------------
 * Solve Poisson's equation for a rectangular box with Dirichlet
 * boundary conditions on each face.
 * source[i, j, k] is accessed with source[((k * ysize) + j) * xsize + i]
 * 
 * @returns:
 *      - uint8_t error: A variable that will be 1 if an error has
 *                       occurred and 0 otherwise.
*/
int poisson_args_t::poissonDirichlet(void)
{
	uint8_t error = EXIT_SUCCESS;

	// Allocate memory for the potenial calculations to be stored
    size_t size = (size_t) y_size * z_size * x_size * sizeof(double); // Calculate the amount of memory needed
	double* input = (double *) malloc(size);

	// Check if memory for 'input' was successfully allocated
	if (!input) {
		fprintf(stderr, "malloc failure\n"); // Print error
		error = ERROR_SYMBOL;

		return error;
	}

	// Perform iterations of Jacobi relaxation
	memcpy(input, source, size); // Copy the source distribution as the input for the first iteration
	for (unsigned int iter = 0; iter < num_iters; iter++) { // Perform iteration of Jacobi relaxation
		for (unsigned int x = 0; x < x_size; x++) { // Iterate through cuboid's x values
			for (unsigned int z = 0; z < z_size; z++) { // Iterate through cuboid's z values
				for (unsigned int y = 0; y < y_size; y++) { // Iterate through cuboid's y values
					double res = 0; // Set the result for this iteration to 0

					if (x < x_size - 1)
						res += input[((z * y_size) + y) * x_size + (x + 1)];
					else
						res += V_bound;
					if (x > 0)
						res += input[((z * y_size) + y) * x_size + (x - 1)];
					else
						res += V_bound;

					if (y < y_size - 1)
						res += input[((z * y_size) + (y + 1)) * x_size + x];
					else
						res += V_bound;
					if (y > 0)
						res += input[((z * y_size) + (y - 1)) * x_size + x];
					else
						res += V_bound;

					if (z < z_size - 1)
						res += input[(((z + 1) * y_size) + y) * x_size + x];
					else
						res += V_bound;
					if (z > 0)
						res += input[(((z - 1) * y_size) + y) * x_size + x];
					else
						res += V_bound; 

					res -= delta * delta * source[((z * y_size) + y) * x_size + x];

					res /= 6;

					potential[((z * y_size) + y) * x_size + x] = res;
				}
			}
		}
		memcpy(input, potential, size); // Copy the calculated potential as the input for the next iteration
	}
	free(input);

	return error;
}
