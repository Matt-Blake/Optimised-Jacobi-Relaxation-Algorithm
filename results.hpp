/* ****************************************************************
 * results.hpp
 *
 * Part of ENCE464 Assignment 2
 * Source file containing the functions needed to save the results
 * of the Algorithim Optimisation program to a .csv file. These
 * results include the time taken to solve Poisson's equation
 * using Jacobi relaxation for different numbers of threads and
 * iterations.
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 03/10/2020
 * ***************************************************************/

#ifndef RESULTS_H
#define RESULTS_H


/*
 * Function:    calculateResults
 * ------------------------------
 * Calculates the time taken for Poisson's equation to be
 * solved using Jacobi relaxation.
 *
 * @params:
 *      - poisson_args_t* poisson_args: The arguments needed to solve
 * 					      				Poisson's equation.
 * ---------------------
 */
void* calculateResults(poisson_args_t* poisson_args);


#endif /* RESULTS_H */