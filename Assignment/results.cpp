/* ****************************************************************
 * main.cpp
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
#include <iostream>
#include <fstream>

#define OUTPUT_FILEPATH		  	  "../Code outputs/Poisson Results.csv"	// The filepath of the .csv to write results to
#define CSV_HEADER_STRING		  "Threads,Time Taken (ms)\n" // A str//ing containing the header text for the .csv results file
#define MAX_STRING_SIZE			  100				// The maximum number of chars to be stored in a string (doesn't include '\0')


/*
 * Function:    saveResultsToCSV
 * ------------------------------
 * Save the results of the Algorithim Optimisation program
 * to a .csv file.
 *
 * @params:
 * 		- 
 * --------------------- 
 */
static void* saveResultsToCSV(void)
{
	std::ofstream output_file;

	// Initalise file
	output_file.open(OUTPUT_FILEPATH);
	output_file << (char*) CSV_HEADER_STRING; // Save table headers to .csv file

	// Calculate the results without using threading
	getNonThreadedSumResults(sum_arguments); 

	// Calculate the results without using threading
	for(int i = 1; i <= sum_arguments->nthreads; i++) { // Iterate through different numbers of threads
		getThreadedSumResults(sum_arguments, i, &output_file);
	}

	output_file.close();

	return NULL;
}
