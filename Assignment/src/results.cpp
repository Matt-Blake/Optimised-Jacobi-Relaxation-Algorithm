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
#include <vector>
#include "poisson.hpp"

#define OUTPUT_FILEPATH		  	  "../../../Code outputs/Assignment/Poisson Results.csv"	// The filepath of the .csv to write results to
#define CSV_HEADER_STRING		  "Size,Iterations,Time Taken (ms)\n" // A str//ing containing the header text for the .csv results file
#define MAX_STRING_SIZE			  100				// The maximum number of chars to be stored in a string (doesn't include '\0')
#define S_TO_NS 			  	  1000000000ULL 	// Conversion factor from seconds to nanoseconds
#define MS_TO_NS				  1000000 			// Conversion factor from milliseconds to nanoseconds


/*
 * Function:    saveTableToCSV
 * ------------------------------
 * Save the results to a .csv file as a table.
 *
 * @params:
 * 		- std::vector<std::string> output_strings: The strings to save.
 * --------------------- 
 */
static void* saveResultsToCSV(std::vector<std::string> output_strings)
{
	std::ofstream output_file;
    int num_results;

	// Initalise file
	output_file.open(OUTPUT_FILEPATH);
	output_file << (char*) CSV_HEADER_STRING; // Save table headers to .csv file

	// Save results
	num_results = output_strings.size(); // Calculate the number of results
    for(int i=0; i < num_results; i++) {
        output_file << output_strings.at(i); // Save result to .csv file
		std::cout << output_strings.at(i);
    }

	output_file.close();

	return NULL;
}


/*
 * Function:    getResultsString
 * ------------------------------
 * Returns a formatted result string based on the time taken to solve
 * Poisson's equation with a particular set of arguments. This
 * string can be printed or saved in a .txt or .csv file.
 *
 * @params:
 *      - poisson_args_t* poisson_args: The arguments needed to solve
 * 					      				Poisson's equation.
 * 		- uint64_t nanoseconds: The time taken to solve Poisson's
 * 								equation (in ns).
 * @returns:
 * 		- std::string: A string containing the results.
 * ---------------------
 */
std::string getResultsString(poisson_args_t* poisson_args, double nanoseconds)
{
	static 
	std::string size_string;
	std::string iterations_string;
	std::string time_string;
	std::string cores_string;
	std::string result_string;

	// Convert results to strings
	size_string = std::to_string(poisson_args->x_size);
	iterations_string = std::to_string(poisson_args->num_iters);
	time_string = std::to_string(nanoseconds);
	cores_string = std::to_string(poisson_args->num_cores);
	

	// Save result
	result_string = size_string + ',' + iterations_string + ',' + cores_string + ',' + time_string + '\n';

	return result_string;

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
    double nanoseconds;
	double milliseconds;
	std::string result_string;
	std::vector<std::string> result_strings;

	// Calculate the time taken to solve Poisson's equation
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
    poisson_dirichlet(poisson_args->source, poisson_args->potential, poisson_args->V_bound, poisson_args->x_size,
					  poisson_args->y_size, poisson_args->z_size, poisson_args->delta, poisson_args->num_iters,
					  poisson_args->num_cores); // Solve Poisson's equation based on arguments
    clock_gettime(CLOCK_MONOTONIC, &end); // End timer
    nanoseconds = (double) (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec); // Calculate time taken
	milliseconds = nanoseconds/MS_TO_NS;
	result_string = getResultsString(poisson_args, milliseconds); // Convert results to formatted string
	result_strings.push_back(result_string); // Save results string in a vector
	saveResultsToCSV(result_strings); // Save all results strings to a .csv file

	return NULL;
}
