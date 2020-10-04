/* ****************************************************************
 * main.cpp
 *
 * Part of ENCE464 Assignment 2
 * Main function of Algorithim Optimisation program. This program
 * attemps to find the potential in a 3-D cube, for a given charge
 * distribution. This is solved through Poisson's equation using
 * Jacobi relaxation. 
 *
 * Based off poisson_test.cpp - Michael Hayes, UC ECE
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 03/10/2020
 * ***************************************************************/


/* TODO LIST for MATT AND DERRICK
- Check if thread working properly
- Suss threading so it works with a non-divisable amount
- Implement loop unrolling
- Implement profiling properly
- Look through lecture notes for other optimisations to implement
- Make sure function is still sussing it properly (check final pointer for memcopy
  and if it should be < z_end or <= z_end)
*/

#include <cstdio>
#include <cstdlib>
#include "poisson.hpp"
#include "results.hpp"


/*
 * Function:    main
 * ------------------------------
 * The main function of Algorithim Optimisation program. This function
 * runs the program which solves Poisson's equation using Jacobi
 * relaxation for a user defined rectangular box.
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
int main(int argc, char* argv[])
{
    poisson_args_t poisson_user_args;
    poisson_args_t* poisson_args_pointer;

    // Calculate the potential in a cuboid based on Poisson's equation
    poisson_args_pointer = &poisson_user_args;
    poisson_args_pointer->initPoissonArgs(argc, argv); // Initalise arguments for solving
    calculateResults(poisson_args_pointer); // Solve Poisson's equation and save the time taken
    
    return EXIT_SUCCESS;
}
