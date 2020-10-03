/* ****************************************************************
 * results.cpp
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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include "poisson.hpp"

#define OUTPUT_FILEPATH		  	  "../Code outputs/Poisson Results.csv"	// The filepath of the .csv to write results to
#define CSV_HEADER_STRING		  "Threads,Time Taken (ms)\n" // A str//ing containing the header text for the .csv results file
#define MAX_STRING_SIZE			  100				// The maximum number of chars to be stored in a string (doesn't include '\0')
#define S_TO_NS 			  	  1000000000ULL 	// Conversion factor from seconds to nanoseconds

/*
 * Function:    saveResultsToCSV
 * ------------------------------
 * Save the results a list of strings to a .csv file.
 *
 * @params:
 * 		- std::string* results_strings: A pointer to the strings
 *                                		to print.
 * --------------------- 
 */
static void* saveResultsToCSV(std::string* output_strings)
{
	std::ofstream output_file;
    int num_results;

	// Initalise file
	output_file.open(OUTPUT_FILEPATH);
	output_file << (char*) CSV_HEADER_STRING; // Save table headers to .csv file

	// Save results
	num_results = sizeof(output_strings); // Calculate the number of results
    for(int i=0; i < num_results; i++) {
        output_file << output_strings[i]; // Save result to .csv file
    }

	output_file.close();

	return NULL;
}


/*
 * Function:    saveResults
 * ------------------------------
 * Saves the time taken for Poisson's equation to be
 * solved using Jacobi relaxation to a .csv.
 *
 * @params:
 *      - poisson_args_t* poisson_args: The arguments needed to solve
 * 					      				Poisson's equation.
 * 		- uint64_t nanoseconds: The time taken to solve Poisson's
 * 								equation (in ns).
 * ---------------------
 */
void* saveResults(poisson_args_t* poisson_args, uint64_t nanoseconds)
{
	std::string size_string;
	std::string iterations_string;
	std::string time_string;

	// Convert results to strings
	size_string = std::to_string(poisson_args->x_size);
	iterations_string = std::to_string(poisson_args->num_iters);
	time_string = std::to_string(nanoseconds);

	// Save result
	std::cout << size_string + ',' + iterations_string + ',' + time_string + '\n';
	//saveResultsToCSV(results);

	return NULL;

}


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
void* calculateResults(poisson_args_t* poisson_args)
{
	struct timespec start, end;
    uint64_t nanoseconds;

	// Calculate the time taken to solve Poisson's equation
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
    poisson_args->poissonDirichlet(); // Solve Poisson's equation based on arguments
    clock_gettime(CLOCK_MONOTONIC, &end); // End timer
    nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec); // Calculate time taken
	saveResults(poisson_args, nanoseconds);

	return NULL;
}