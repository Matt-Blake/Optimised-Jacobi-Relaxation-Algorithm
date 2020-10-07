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
#include <thread>
#include <vector>
#include "poisson.hpp"

#include <iostream>

#define VOXEL_SPACING       0.1     // The spacing (in all directions) between voxels (meters)
#define V_BOUND             0       // Voltage potential on the box boundary
#define CHARGE_VALUE        1       // The value of a point charge in the charge distribution
#define ERROR_SYMBOL        1       // The value to return if an error occurs
#define THREADS_PER_CORE	1		// The number of threads to use per core (optimal performance at 1)


/*
 * Class:    thread_args
 * ------------------------------
 * Contains the information needed to split up Jacobi
 * relaxation into slices solvable by threads. This is
 * done by splitting the z-axis of the cuboid into slices.
 * 
 * 
 * @members:
        - unsigned int thread_number: A number corresponding to the thread
 * ---------------------
 */
typedef class thread_args
{
    public:
        unsigned int thread_number;
} thread_args_t;


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
 * Perform an iteration of Jacobi Relaxation using a single thread.
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
 * 		- unsigned int x_size: The number of x-axis elements.
 * 		- unsigned int y_size: The number of y-axis elements.
 * 		- unsigned int z_size: The number of z-axis elements.
 *		- double delta: The spacing between voxels.
 *		- unsigned int num_threads: The number of threads being used.
 *		- thread_args_t* thread: A pointer to the start and finishing z-axis
 *								  values to iterate through for each thread
 * --------------------- 
*/
void* performJacobiIteration(double* input, double* potential, double* source, 
							   unsigned int x_size, unsigned int y_size, unsigned int z_size,
							   double delta, unsigned int num_threads, thread_args_t* thread)
{
	// Iterate through each voxel, calculating the potential via Jacobi's relaxation
	for (unsigned int z = thread->thread_number; z < (z_size - 1); z += num_threads) { // Iterate through cuboid's z values
		for (unsigned int y = 1; y < (y_size - 1); y++) { // Iterate through cuboid's y values
			for (unsigned int x = 1; x < (x_size - 1); x++) { // Iterate through cuboid's x values
				
				double result = 0;
				
				result += input[((z * y_size) + y) * x_size + (x + 1)]; // V[x+1, y, z]
				result += input[((z * y_size) + y) * x_size + (x - 1)]; // V[x-1, y, z]
				result += input[((z * y_size) + (y + 1)) * x_size + x]; // V[x, y+1, z]
				result += input[((z * y_size) + (y - 1)) * x_size + x]; // V[x, y-1, z]
				result += input[(((z + 1) * y_size) + y) * x_size + x]; // V[x, y, z+1]
				result += input[(((z - 1) * y_size) + y) * x_size + x]; // V[x, y, z-1]
				result -= delta * delta * source[((z * y_size) + y) * x_size + x]; // Subtract the effect of the source
				result /= 6;

				potential[((z * y_size) + y) * x_size + x] = result; // Store potential result for current voxel
			}
		}
	}
	
	return NULL;
}


/*
 * Function:    performJacobiBoundaryIteration
 * ------------------------------
 * Perform an iteration of Jacobi Relaxation by focusing on
 * boundary conditions.
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
 *		- unsigned int num_threads: The number of threads being used.
 *		- thread_args_t* thread: A pointer to the start and finishing axis values
 *								  to iterate through for each thread
 * --------------------- 
*/
void* performJacobiBoundaryIteration(double* input, double* potential, double* source, double V_bound,
									unsigned int x_size, unsigned int y_size, unsigned int z_size,
									double delta, unsigned int num_threads, thread_args_t* thread)
{
	unsigned int x_bound = 0;
	unsigned int y_bound = 0;
	unsigned int z_bound = 0;

	// Calculate x-axis boundary cases
	for (unsigned int z = thread->thread_number; z < z_size; z += num_threads) { // Iterate through cuboid's z values
		for (unsigned int y = 0; y < y_size; y++) { // Iterate through cuboid's y values
			
			// x-axis 0 boundary case
			double result = 0;
			x_bound = 0;
			result += V_bound; // V[x-1, y, z]
			result += input[((z * y_size) + y) * x_size + (x_bound + 1)]; // V[x+1, y, z]
			if (y < (y_size - 1)) {
				result += input[((z * y_size) + (y + 1)) * x_size + x_bound]; // V[x, y+1, z]
			} else {
				result += V_bound;
			}
			if (y > 0) {
				result += input[((z * y_size) + (y - 1)) * x_size + x_bound]; // V[x, y-1, z]
			} else {
				result += V_bound;
			}
			if (z < (z_size - 1)) {
				result += input[(((z + 1) * y_size) + y) * x_size + x_bound]; // V[x, y, z+1]
			} else {
				result += V_bound;
			}
			if (z > 0) {
				result += input[(((z - 1) * y_size) + y) * x_size + x_bound]; // V[x, y, z-1]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z * y_size) + y) * x_size + x_bound]; // Subtract the effect of the source
			result /= 6;
			potential[((z * y_size) + y) * x_size + x_bound] = result; // Store potential result for current voxel
			
			// x-axis max boundary case
			result = 0;
			x_bound = x_size - 1;
			result += V_bound; // V[x+1, y, z]
			result += input[((z * y_size) + y) * x_size + (x_bound - 1)]; // V[x-1, y, z]
			if (y < (y_size - 1)) {
				result += input[((z * y_size) + (y + 1)) * x_size + x_bound]; // V[x, y+1, z]
			} else {
				result += V_bound;
			}
			if (y > 0) {
				result += input[((z * y_size) + (y - 1)) * x_size + x_bound]; // V[x, y-1, z]
			} else {
				result += V_bound;
			}
			if (z < (z_size - 1)) {
				result += input[(((z + 1) * y_size) + y) * x_size + x_bound]; // V[x, y, z+1]
			} else {
				result += V_bound;
			}
			if (z > 0) {
				result += input[(((z - 1) * y_size) + y) * x_size + x_bound]; // V[x, y, z-1]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z * y_size) + y) * x_size + x_bound]; // Subtract the effect of the source
			result /= 6;
			potential[((z * y_size) + y) * x_size + x_bound] = result; // Store potential result for current voxel
		}
	}

	// Calculate y-axis boundary cases
	for (unsigned int z = thread->thread_number; z < z_size; z += num_threads) { // Iterate through cuboid's z values
		for (unsigned int x = 0; x < x_size; x++) { // Iterate through cuboid's x values
			
			// y-axis 0 boundary case
			double result = 0;
			y_bound = 0;
			result += V_bound; // V[x, y-1, z]
			result += input[((z * y_size) + (y_bound + 1)) * x_size + x]; // V[x, y+1, z] 
			if (x < (x_size - 1)) {
				result += input[((z * y_size) + y_bound) * x_size + (x + 1)]; // V[x+1, y, z]
			} else {
				result += V_bound;
			}
			if (x > 0) {
				result += input[((z * y_size) + y_bound) * x_size + (x - 1)]; // V[x-1, y, z]
			} else {
				result += V_bound;
			}
			if (z < (z_size - 1)) {
				result += input[(((z + 1) * y_size) + y_bound) * x_size + x]; // V[x, y, z+1]
			} else {
				result += V_bound;
			}
			if (z > 0) {
				result += input[(((z - 1) * y_size) + y_bound) * x_size + x]; // V[x, y, z-1]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z * y_size) + y_bound) * x_size + x]; // Subtract the effect of the source
			result /= 6;
			potential[((z * y_size) + y_bound) * x_size + x] = result; // Store potential result for current voxel

			// y-axis max boundary case
			result = 0;
			y_bound = y_size - 1;
			result += V_bound; // V[x, y+1, z]
			result += input[((z * y_size) + (y_bound - 1)) * x_size + x]; // V[x, y-1, z] 
			if (x < (x_size - 1)) {
				result += input[((z * y_size) + y_bound) * x_size + (x + 1)]; // V[x+1, y, z]
			} else {
				result += V_bound;
			}
			if (x > 0) {
				result += input[((z * y_size) + y_bound) * x_size + (x - 1)]; // V[x-1, y, z]
			} else {
				result += V_bound;
			}
			if (z < (z_size - 1)) {
				result += input[(((z + 1) * y_size) + y_bound) * x_size + x]; // V[x, y, z+1]
			} else {
				result += V_bound;
			}
			if (z > 0) {
				result += input[(((z - 1) * y_size) + y_bound) * x_size + x]; // V[x, y, z-1]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z * y_size) + y_bound) * x_size + x]; // Subtract the effect of the source
			result /= 6;
			potential[((z * y_size) + y_bound) * x_size + x] = result; // Store potential result for current voxel
		}
	}

	// Calculate z-axis boundary cases
	for (unsigned int y = thread->thread_number; y < y_size; y += num_threads) { // Iterate through cuboid's y values
		for (unsigned int x = 0; x < x_size; x++) { // Iterate through cuboid's x values
			
			// z-axis 0 boundary case
			double result = 0;
			z_bound = 0;
			result += V_bound; // V[x, y, z-1]
			result += input[(((z_bound + 1) * y_size) + y) * x_size + x];; // V[x, y, z+1]
			if (x < (x_size - 1)) {
				result += input[((z_bound * y_size) + y) * x_size + (x + 1)]; // V[x+1, y, z]
			} else {
				result += V_bound;
			}
			if (x > 0) {
				result += input[((z_bound * y_size) + y) * x_size + (x - 1)]; // V[x-1, y, z]
			} else {
				result += V_bound;
			}
			if (y < (y_size - 1)) {
				result += input[(((z_bound) * y_size) + (y + 1)) * x_size + x]; // V[x, y+1, z]
			} else {
				result += V_bound;
			}
			if (y > 0) {
				result += input[(((z_bound) * y_size) + (y - 1)) * x_size + x]; // V[x, y-1, z]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z_bound * y_size) + y) * x_size + x]; // Subtract the effect of the source
			result /= 6;
			potential[((z_bound * y_size) + y) * x_size + x] = result; // Store potential result for current voxel

			// z-axis max boundary case
			result = 0;
			z_bound = z_size - 1;
			result += V_bound; // V[x, y, z+1]
			result += input[(((z_bound - 1) * y_size) + y) * x_size + x];; // V[x, y, z-1]
			if (x < (x_size - 1)) {
				result += input[((z_bound * y_size) + y) * x_size + (x + 1)]; // V[x+1, y, z]
			} else {
				result += V_bound;
			}
			if (x > 0) {
				result += input[((z_bound * y_size) + y) * x_size + (x - 1)]; // V[x-1, y, z]
			} else {
				result += V_bound;
			}
			if (y < (y_size - 1)) {
				result += input[((z_bound * y_size) + (y + 1)) * x_size + x]; // V[x, y+1, z]
			} else {
				result += V_bound;
			}
			if (y > 0) {
				result += input[((z_bound * y_size) + (y - 1)) * x_size + x]; // V[x, y-1, z]
			} else {
				result += V_bound;
			}
			result -= delta * delta * source[((z_bound * y_size) + y) * x_size + x]; // Subtract the effect of the source
			result /= 6;
			potential[((z_bound * y_size) + y) * x_size + x] = result; // Store potential result for current voxel
		}
	}

	return NULL;
}

/*
 * Function:    threadJacobiIteration
 * ------------------------------
 * Perform an iteration of Jacobi Relaxation. This is done by splitting
 * the z-axis of the cubiod into threads.
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
 *		- unsigned int num_threads: The number of threads being used.
 * --------------------- 
*/
void* threadJacobiIteration(double* input, double* potential, double* source, double V_bound, 
						   unsigned int x_size, unsigned int y_size, unsigned int z_size,
						   double delta, unsigned int num_threads)
{
	thread_args_t threads[num_threads]; // Create an array of thread argument structs
	std::vector<std::thread> thread_vector(num_threads);
	thread_args_t threads_boundary[num_threads]; // Create an array of thread argument structs
	std::vector<std::thread> thread_boundary_vector(num_threads);

	// Spawn thread to do a sub-section of Jacobi relaxation
	for (int i = 0; i < num_threads; i++) {
		threads[i].thread_number = i + 1;
		thread_vector[i] = std::thread(performJacobiIteration, input, potential, source, x_size, y_size,
								       z_size, delta, num_threads, &threads[i]); // Add thread to vector
	}
	
	// Resolve threads once they are finished
	for (int i = 0; i < num_threads; i++) {
		thread_vector[i].join();
	}

	// Spawn threads to do the boundary cases of Jacobi relaxation
	for (int i = 0; i < num_threads; i++) {
		threads_boundary[i].thread_number = i;
		thread_boundary_vector[i] = std::thread(performJacobiBoundaryIteration, input, potential, source,
												V_bound, x_size, y_size, z_size, delta, num_threads,
												&threads_boundary[i]); // Add thread to vector
	}

	// Resolve threads once they are finished
	for (int i = 0; i < num_threads; i++) {
		thread_boundary_vector[i].join();
	}

	return NULL;
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
	// Allocate memory for the potential calculations to be stored
    size_t size = (size_t) xsize * ysize * zsize * sizeof(double); // Calculate the amount of memory needed
	double* input = (double*) malloc(size);
	double* temp;
	unsigned int num_threads;
	
	// Check if memory for 'input' was successfully allocated
	if (!input) {
		fprintf(stderr, "malloc failure\n"); // Print error

		return;
	}

	// Calculate number of threads needed to split up Jacobi relaxation into slices
	num_threads = numcores * THREADS_PER_CORE;
	if (num_threads > zsize) { // Catch case where the cubiod is very small, so data can't be split into threads
		num_threads = zsize;
	}

	// Perform iterations of Jacobi relaxation
	memcpy(input, source, size); // Copy the source distribution as the input for the first iteration
	for (unsigned int iter = 0; iter < numiters; iter++) {
		
		// Split up the z-axis values and spawn the threads to do an iteration of Jacobi relaxation
		threadJacobiIteration(input, potential, source, Vbound, xsize, ysize, zsize, delta, num_threads);

		// Swap pointers to prepare fo the next iteration
		if (iter != (numiters - 1)) { // If another iteration will occur
			temp = input;
			input = potential; // Copy the calculated potential as the input for the next iteration
			potential = temp; // Copy the potential that will be modified to a different address as the input
		}
	}
	
	// For checking potential is calculated correctly (DON'T DELETE YET, i'll delete it at the end)
	for (int i=0; i < zsize; i++) {
		for (int j=0; j < ysize; j++) {
			for (int k=0; k < xsize; k++) {
				if (k == 0) {
					std::cout << '[';
				} 
				std::cout << potential[((k * ysize) + j) * xsize + i];
				
				if (k == (xsize - 1)) {
					std::cout << ']';
				} else {
					std::cout << ", ";
				}
				//std::cout << "potential ";
				//std::cout << potential[((k * ysize) + j) * xsize + i];
				//std::cout << '\n';
				//std::cout << "z value ";
				//std::cout << i;
	    		//std::cout << '\n';
				//std::cout << "y value ";
				//std::cout << j;
	    		//std::cout << '\n';
				//std::cout << "x value ";
				//std::cout << k;
	    		//std::cout << '\n';
				//std::cout << '\n';
				
			}
		}
	}
	
	free(input);
}
