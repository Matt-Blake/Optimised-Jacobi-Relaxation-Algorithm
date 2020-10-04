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

#include <cstdio>
#include <cstdlib>
#include <cstring>
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
 * --------------------- 
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
 * --------------------- 
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
 * --------------------- 
 */
void* poisson_args_t::initPoissonArgs(int argc, char** argv)
{  
    allocateUserInputs(argc, argv); // Initalise members based on user inputs
    allocateVolume(); // Dynamically create space for the 3-D volume arrays

    return NULL;
}


/*
 * Function:    performJacobiIteration
 * ------------------------------
 * Perform an iteration of Jacobi Relaxation
 * input[i, j, k] is accessed with input[((k * y_size) + j) * x_size + i].
 * This function is laid out to fulfill the ENCE464 Assignment 2
 * specifications. However, future program develpments would have this
 * function as a method to the poisson_args_t object.
 * 
 * @params:
 *      - double* input: A 3-D array of potential values representing the
 * 						 previous iteration  of calculations using
 * 						 Jacobi relaxation.
 * 		- double* potential: A 3-D array of potential values representing the
 * 						    the calculations that will be performed this
 * 						    iteration of Jacobi relaxation.
 * 		- double* source: A 3-D array representing the charge distribution.
 * 		- double V-bound: The potential on the boundary conditions.
 * 		- unsigned int x_size: The number of x-axis elements.
 * 		- unsigned int y_size: The number of y-axis elements.
 * 		- unsigned int z_size: The number of z-axis elements.
 *		- double delta: The spacing between voxels.
 * @returns:
 *      - double* potential: A 3-D array of potential values representing the
 * 						    the calculations performed this iteration of
 * 						    Jacobi relaxation.
 * --------------------- 
*/
double* performJacobiIteration(double* input, double* potential, double* source, double V_bound, 
							unsigned int x_size, unsigned int y_size, unsigned int z_size, double delta)
{
	// Iterate through each voxel, calculating the potential via Jacobi's relaxation
	for (unsigned int z = 0; z < z_size; z++) { // Iterate through cuboid's x values
		for (unsigned int y = 0; y < y_size; y++) { // Iterate through cuboid's z values
			for (unsigned int x = 0; x < x_size; x++) { // Iterate through cuboid's y values
				
				double res = 0; // Initalise the result for the current voxel to 0

				// Calculate V[x+1, y, z, iter]
				if (x < x_size - 1)
					res += input[((z * y_size) + y) * x_size + (x + 1)];
					
				else // Boundary condition
					res += V_bound;

				// Calculate V[x-1, y, z, iter]
				if (x > 0)
					res += input[((z * y_size) + y) * x_size + (x - 1)];
					
				else // Boundary condition
					res += V_bound;

				// Calculate V[x, y+1, z, iter]
				if (y < y_size - 1)
					res += input[((z * y_size) + (y + 1)) * x_size + x];
					
				else // Boundary condition
					res += V_bound;
				
				// Calculate V[x, y-1, z, iter]
				if (y > 0)
					res += input[((z * y_size) + (y - 1)) * x_size + x];
					
				else // Boundary condition
					res += V_bound;

				// Calculate V[x, y, z+1, iter]
				if (z < z_size - 1)
					res += input[(((z + 1) * y_size) + y) * x_size + x];
					
				else // Boundary condition
					res += V_bound;

				// Calculate V[x, y, z-1, iter]
				if (z > 0)
					res += input[(((z - 1) * y_size) + y) * x_size + x];
				else // Boundary condition
					res += V_bound; 

				// Subtract the effect of the source
				res -= delta * delta * source[((z * y_size) + y) * x_size + x];

				// Divide result by 6 as per Jacobi's relaxation
				res /= 6;
			
				// Store potential result for current voxel
				potential[((z * y_size) + y) * x_size + x] = res;
			}
		}
	}

	return potential;
}


/*
 * Function:    poisson_dirichlet
 * ------------------------------
 * Solve Poisson's equation for a rectangular box with Dirichlet
 * boundary conditions on each face. The remaining documentation
 * is included to fulfill the ENCE464 Assignment 2 specifications.
 * However, future program develpments would have this function
 * as a method to the poisson_args_t object.
 * 
 * --------------------- 
*/
/// Solve Poisson's equation for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param Vbound is the potential on the boundary
/// \param xsize is the number of elements in the x-direction
/// \param ysize is the number of elements in the y-direction
/// \param zsize is the number of elements in the z-direction
/// \param delta is the voxel spacing in all directions
/// \param numiters is the number of iterations to perform
/// \param numcores is the number of CPU cores to use.  If 0, an optimal number is chosen
void poisson_dirichlet (double * __restrict__ source,
                        double * __restrict__ potential,
                        double Vbound,
                        unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta,
                        unsigned int numiters, unsigned int numcores)

{
	uint8_t error = EXIT_SUCCESS;

	// Allocate memory for the potential calculations to be stored
    size_t size = (size_t) ysize * zsize * xsize * sizeof(double); // Calculate the amount of memory needed
	double* input = (double*) malloc(size);
	double* temp;

	// Check if memory for 'input' was successfully allocated
	if (!input) {
		fprintf(stderr, "malloc failure\n"); // Print error
		error = ERROR_SYMBOL;

		return;
	}

	// Perform iterations of Jacobi relaxation
	input = source; // Copy the source distribution as the input for the first iteration
	for (unsigned int iter = 0; iter < numiters; iter++) {
		potential = performJacobiIteration(input, potential, source, Vbound, xsize, ysize,
										   zsize, delta); // Perform iteration of Jacobi relaxation
		
		// Swap pointers to prepare fo the next iteration
		if (iter != (numiters - 1)) { // If another iteration will occur
			temp = input;
			input = potential; // Copy the calculated potential as the input for the next iteration
			potential = temp; // Copy the potential that will be modified to a different address as the input
		}
	}

	free(input);
}
