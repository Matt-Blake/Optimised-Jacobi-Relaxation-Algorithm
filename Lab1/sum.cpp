/* ****************************************************************
 * sum.cpp
 *
 * Part of ENCE464 Assignment 2
 * This program uses a specified number of threads to sum up
 * a specified number of random values.
 *
 * Based off sum.c - Andre Renaud, UC ECE Department, 2020
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 27/09/2020
 *
 * ***************************************************************/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
#include <thread>
#include <iostream>
#include <fstream>   

#define S_TO_NS 			  	  1000000000ULL 	// Conversion factor from seconds to nanoseconds
#define MS_TO_NS 			  	  1000000ULL		// Conversion factor from millseconds to nanoseconds
#define EXIT_STATUS			  	  1					// The status value returned when the program is exited due to a malloc error
#define COUNT_ARG_POS 		  	  1					// The position of the count argument to main
#define THREADS_ARG_POS 	  	  2					// The position of the thread argument to main

#define BASE_VALUES_TO_SUM 	  	  10000 			// The number of values to sum if no inputs are given to main
#define BASE_NUM_THREADS 	  	  10   				// The number of threads to use if no inputs are given to main
#define MAX_VALUE			  	  100				// The values will being summed will be between 0 and (MAX_VALUE - 1)
#define OVERWRITE_VALUES		  true				// If the input number of summed values should be overwritten and divisable by the number of threads

#define OUTPUT_FILEPATH		  	  "../Code outputs/Thread Results.csv"	// The filepath of the .csv to write results to
#define CSV_HEADER_STRING		  "Threads,Time Taken (ms)\n" // A str//ing containing the header text for the .csv results file
#define MAX_STRING_SIZE			  100				// The maximum number of chars to be stored in a string (doesn't include '\0')


/*
 * Class:    sum_args
 * ------------------------------
 * Defines stucture to contain the information needed
 * for summation results.
 * 
 * @members:
 *  	- int nthreads: The maximum number of threads to use
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 * @methods:
 * 		- allocateValues: Fills 'vals' with random values
 * 		- calcNumValues: Calculates the number of values to
 * 						 be summed
 * ---------------------
 */
class sum_args
{
	public:
		double* vals;			  
		size_t count;
		double nthreads;
		void* allocateValues(void);
	private:
		void* calcNumValues(void);
};


/*
 * Function:    calcNumValues
 * ------------------------------
 * Calculates the minimum number of values that allow
 * the summed values can be divided equally by all of
 * the numbers of threads being compared.
 *
 * @updates:
 *      - size_t count: The number of values to be created
 * --------------------- 
 */
void* sum_args::calcNumValues(void)
{
	// Calculate a count that is a multiple of all thread numbers
	count = 1;
	for(int i = 1; i <= nthreads; i++) {
		count *= i;
	}

	return NULL;
}


/*
 * Function:    allocateValues
 * ------------------------------
 * Dynamically allocates 'count' random values to an array.
 *
 * @params:
 *      - int nthreads: The maximum number of threads to use
 * @updates:
 *      - double* vals: A pointer to the allocated values
 * --------------------- 
 */
void* sum_args::allocateValues(void)
{
	double* values;
	
	// If the user inputted value needs to be divsabled by nthreads
	if (OVERWRITE_VALUES) {
		calcNumValues();
	}
	
	// Dynamically allocate space for values
	values = (double*) malloc(sizeof(double) * count); 
	if (!values) // Exit if no space could be dynamically allocated
		exit(1);

	// Assign a random value to each value
	for (size_t i = 0; i < count; i++) { 
		values[i] = rand() % MAX_VALUE;
	}
	vals = values;
	
	return NULL;
}


/*
 * Structure:    thread_args
 * ------------------------------
 * Defines stucture to contain the information needed
 * for arguments to the sum.cpp thread functions.
 * 
 * @members:
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 *      - double thread_sum: The sum of values in 'vals' 
 *      - std::thread* thread: A pointer to the thread used
 * ---------------------
 */
struct thread_args
{
	double* vals;			  
	size_t count;			   
	double thread_sum;
	std::thread* thread;
};


/*
 * Function:    noThreadSum
 * ------------------------------
 * no_thread_sum returns the sum of the values in 'vals', of length 'count'
 * without using threading. This provides a basis for comparison with 
 * threaded methods.
 *
 * @params:
 *      - sum_args sum_arguments: Contains the information
 * 								  needs for the current sum.
 * @return:
 *      - double sum: The sum of the values in 'vals'
 * ---------------------
 */
static double noThreadSum(sum_args* sum_arguments)
{
	double sum = 0;

	// Iterate through the values in 'vals' summing them
	for (size_t i = 0; i < sum_arguments->count; i++) {
		sum += sum_arguments->vals[i];
	}

	return sum;
}


/*
 * Function:    getNonThreadedSumResults
 * ------------------------------
 * Sums the value of an array without using threading. This 
 * sum is then printed.
 *
 * @params:
 * 		- sum_args sum_arguments: Contains the information
 * 								  needs for the current sum.
 * --------------------- 
 */
static void* getNonThreadedSumResults(sum_args* sum_arguments)
{
	struct timespec start, end;
	uint64_t nanoseconds;
	double sum;

	// Calculate the sum of 'vals' and the time taken to perform this sum without using threads
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
	sum = noThreadSum(sum_arguments);
	clock_gettime(CLOCK_MONOTONIC, &end); // End timer
	nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec);

	// Print results
	printf("no thread sum: %f\nTook %" PRIu64 " ms for %zu iterations\n", sum, nanoseconds / MS_TO_NS, sum_arguments->count); // Print to terminal

	return NULL;
}


/*
 * Function:    threadSum
 * ------------------------------
 * thread_sum_func uses a thread to sum the the sum of the values. 
 *
 * @params:
 *      - thread_args* thread_arguments: A pointer to the arguments needed
 *  								     to initalise a thread_args struct
 * ---------------------
 */
static void threadSum(void)//thread_args* thread_arguments)
{
	//double sum = 0;

	// Iterate through the values in 'thread_arg->vals' summing them
	//for (size_t i = 0; i < thread_arguments->count; i++) {
	//	sum += thread_arguments->vals[i];
	//}
	//thread_arguments->thread_sum = sum; // Store the sum in the thread structure
}


/*
 * Function:    sumValues
 * ------------------------------
 * Break up 'vals' into a series of slices and spawn a thread to sum each
 * slice, then return the sum of these intermediate sums. For simplicity,
 * each thread sums an equal number of values.
 *
 * @params:
 *      - sum_args sum_arguments: Contains the information
 * 								  needs for the current sum.
 *      - int nthreads: The number of threads to use
 * @return:
 *      - double sum: The sum of values
 * --------------------- 
 */
static double sumValues(sum_args* sum_arguments, int nthreads)
{
	thread_args thread_arg[nthreads];
	double sum = 0;

	// Print error and exit the number of values can not be split equally amongst the threads
	if (sum_arguments->count % nthreads) {
		fprintf(stderr, "count %zu must be divisible by %d\n", sum_arguments->count, nthreads);
		return 0;
	}

	// Split up the values in 'vals' equally and iteratively spawn the threads
	for (int i = 0; i < nthreads; i++) {

		// Initalise thread_arg structure
		thread_arg[i].count = (sum_arguments->count) / nthreads;
		thread_arg[i].vals = &(sum_arguments->vals[i * (sum_arguments->count) / nthreads]);
		thread_arg[i].thread_sum = 0;
		std::thread summing_thread(threadSum);//, &thread_arg[i]); // ERROR OCCOURING HERE
		//thread_arg[i].thread = &summing_thread;
	}

	// Wait for each thread to finish, then add in its partial sum to the total sum
	for (int i = 0; i < nthreads; i++) {
		//thread_arg[i].thread->join();
		//delete thread_arg[i].thread;
		//sum += thread_arg[i].thread_sum;
	}


	return sum;
}


/*
 * Function:    getThreadedSumResults
 * ------------------------------
 * Sums the value of an array using threading. The results
 * are then saved to a .csv file and printed to the terminal.
 *
 * @params:
 * 		- sum_args sum_arguments: Contains the information
 * 								  needs for the current sum.
 *      - int nthreads: The number of threads to use
 * 		- std::ofstream* output_file: A pointer to the opened output file
 * --------------------- 
 */
static char* getThreadedSumResults(sum_args* sum_arguments, int nthreads, std::ofstream* output_file)
{
	struct timespec start, end;
	uint64_t nanoseconds;
	double sum;
	char* results_string = (char*) malloc((MAX_STRING_SIZE + 1) * sizeof(char));
	
	// Calculate the sum of 'vals' and the time taken to perform this sum using threads
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
	sum = sumValues(sum_arguments, nthreads);
	clock_gettime(CLOCK_MONOTONIC, &end); // End timer
	nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec);

	// Print the results
	sprintf(results_string, "%d,%llu\n", nthreads, nanoseconds/MS_TO_NS);  // Save .csv results in a string
	*output_file << results_string; // Write results to .csv file
	printf("thread sum: %f (%d threads)\nTook %" PRIu64 " ms for %zu iterations\n",
	        sum, nthreads, nanoseconds / MS_TO_NS, sum_arguments->count); // Print to terminal

	return NULL;
}


/*
 * Function:    getResults
 * ------------------------------
 * Calculate the time taken to sum summation's arguments 
 * random values. These results are calculated for all number
 * of threads up to 'nthreads' and without using threads. 
 * These results are then printed and saved to a .csv file.
 *
 * @params:
 * 		- sum_args sum_arguments: Contains the information
 * 								  needs for the current sum.
 * --------------------- 
 */
static void* getResults(sum_args* sum_arguments)
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


/*
 * Function:    main
 * ------------------------------
 * The main function of sum.cpp. This function sums a user defined
 * number of random values using a user defined number of threads.
 * This is timed and compared to a case which does not use threading.
 *
 * @params:
 *      - int argc: The number of arguments passed into main
 *      - char** argv: The vector of arguments passed into main.
 * 					   This should contain the number of values
 * 					   to sum, followed by the number of threads
 * 					   to use.
 * ---------------------
 */
int main(int argc, char** argv)
{
	// Initalise variables
	size_t count = BASE_VALUES_TO_SUM;
	int nthreads = BASE_NUM_THREADS;
	double* vals;
	sum_args sum_arguments;

	// Re-initalise variables based on the program's arguments
	if (argc > COUNT_ARG_POS) {
		sum_arguments.count = atoll(argv[COUNT_ARG_POS]);
	}
	if (argc > THREADS_ARG_POS) {
		sum_arguments.nthreads = atoi(argv[THREADS_ARG_POS]);
	}
	
	// Calculate results
	sum_arguments.allocateValues(); // Allocate the values to sum
	getResults(&sum_arguments); // Sum values and store the results

	free(vals); 

	return 0;
}
