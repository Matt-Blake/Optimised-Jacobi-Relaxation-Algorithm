/* ****************************************************************
 * poisson_test.cpp
 *
 * Part of ENCE464 Assignment 2
 * This program attemps to find the potential in a 3-D cube, for a
 * given charge distribution. This is solved through Poisson's
 * equation using Jacobi relaxation. 
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
#define V_BOUND             0       // Voltage potential on the box boundary
#define CHARGE_VALUE        1       // The value of a point charge in the charge distribution
#define ERROR_SYMBOL        1       // The error value to return if not enough inputs were provided

/*
 * Class:    poisson_args
 * ------------------------------
 * Contains the information needed to solve Poisson's
 * equation.
 * 
 * @members:
        - unsigned int num_iters: Number of iterations of Jacobi relaxation to perform
        - unsigned int num_cores: Number of CPU cores to use
        - unsigned int x_size: Size of the rectangular box used in the x dimension
        - unsigned int y_size: Size of the rectangular box used in the y dimension
        - unsigned int z_size: Size of the rectangular box used in the z dimension    
        - double delta: The spacing (in all directions) between voxels (meters)
        - double V_bound: The boundary condition potential
    @methods:
        - initPoissonArgs() - Initialises a poisson_args object
        - allocateUserInputs() - Allocates the user defined members during initalisation
        - allocateVolume() - Allocates memory for the 3-D cuboid's volume
 * ---------------------
 */
typedef class poisson_args
{
	public:
        unsigned int num_iters;
        unsigned int num_cores;    
        unsigned int x_size;
        unsigned int y_size;
        unsigned int z_size;    
        double delta = VOXEL_SPACING;
        double V_bound = V_BOUND;
        double* source;
        double* potential;
        void* initPoissonArgs(int, char**);
    private:
        int allocateUserInputs(int, char**);
        void* allocateVolume(void);
} poisson_args_t;


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
 */
int poisson_args_t::allocateUserInputs(int argc, char** argv)
{  
    unsigned int N;

    // Check if enough user inputs were given
    if (argc < 3)
    {
        fprintf (stderr, "Usage: %s size numiters\n", argv[0]); // Print error
        
        return ERROR_SYMBOL;
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
 * This object can then be passed into the poisson_dirichlet
 * function for solving.
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
 * Function:    main
 * ------------------------------
 * The main function of poisson_test.cpp. This function solves Poisson's
 * equation using Jacobi relaxation for a user defined rectangular box.
 *
 * @params:
 *      - int argc: The number of arguments passed into main
 *      - char** argv: The vector of arguments passed into main.
 * 					   This should contain the 3-D size of the box,
 *                     followed by the number of Jacobi relaxiation
 *                     iterations to perform, then an optional value
 *                     specifying the number of CPU cores to use.         
 * ---------------------
 */
int main (int argc, char* argv[])
{
    poisson_args_t poisson_user_args;
    poisson_args_t* poisson_args_pointer;

    // Calculate the potential in a cuboid based on Poisson's equation
    poisson_args_pointer = &poisson_user_args;
    poisson_args_pointer->initPoissonArgs(argc, argv); // Initalise arguments for solving
    //poisson_dirichlet(poisson_user_args); // Solve Poisson's equation based on arguments
    
    return 0;
}
