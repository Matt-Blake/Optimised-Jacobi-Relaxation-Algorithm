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
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
//#include <thread>
#include <pthread.h>
#include <iostream>
#include <fstream>   

#define S_TO_NS 			  	  1000000000ULL 	// Conversion factor from seconds to nanoseconds
#define MS_TO_NS 			  	  1000000ULL		// Conversion factor from millseconds to nanoseconds
#define BASE_VALUES_TO_SUM 	  	  10000 			// The number of values to sum if no inputs are given to main
#define BASE_NUM_THREADS 	  	  10   				// The number of threads to use if no inputs are given to main
#define COUNT_ARG_POS 		  	  1					// The position of the count argument to main
#define THREADS_ARG_POS 	  	  2					// The position of the thread argument to main
#define MAX_VALUE			  	  100				// The values will being summed will be between 0 and (MAX_VALUE - 1)
#define EXIT_STATUS			  	  1					// The status value returned when the program is exited due to a malloc error

#define OUTPUT_FILEPATH		  	  "../Code outputs/Thread Results.csv"	// The filepath of the .csv to write results to
#define MAX_STRING_SIZE			  100				// The maximum number of chars to be stored in a string (doesn't include '\0')
#define MAX_STRINGS_TO_SAVE   	  100				// The maximum number of strings to be save to 'OUTPUT_FILEPATH' (doesn't include '\0')
#define THREADED_RESULT_ORDER  	  0					// The order in the .csv where the threaded results string should be saved
#define NON_THREADED_RESULT_ORDER 1					// The order where the .csv where the non-threaded results string should be saved

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
 *      - std::thread::id thread_ID: The thread ID used
 * ---------------------
 */
struct thread_args
{
	double* vals;			  
	size_t count;			   
	double thread_sum;
	//std::thread::id thread_ID;
	pthread_t thread;
};


/*
 * Function:    noThreadSum
 * ------------------------------
 * no_thread_sum returns the sum of the values in 'vals', of length 'count'
 * without using threading. This provides a basis for comparison with 
 * threaded methods.
 *
 * @params:
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 * @return:
 *      - double sum: The sum of the values in 'vals'
 * ---------------------
 */
static double noThreadSum(double* vals, size_t count)
{
	double sum = 0;

	// Iterate through the values in 'vals' summing them
	for (size_t i = 0; i < count; i++) {
		sum += vals[i];
	}

	return sum;
}


/*
 * Function:    threadSum
 * ------------------------------
 * thread_sum_func uses a thread to sum the the sum of the values. 
 *
 * @params:
 *      - void* args: A pointer to the arguments needed to initalise
 * 					  a thread_args struct
 * ---------------------
 */
static void* threadSum(void* args)
{
	struct thread_args* thread_arg = (struct thread_args*) args;
	double sum = 0;

	// Iterate through the values in 'thread_arg->vals' summing them
	for (size_t i = 0; i < thread_arg->count; i++) {
		sum += thread_arg->vals[i];
	}
	thread_arg->thread_sum = sum; // Store the sum in the thread structure

	return NULL;
}


/*
 * Function:    sumValues
 * ------------------------------
 * Break up 'vals' into a series of slices and spawn a thread to sum each
 * slice, then return the sum of these intermediate sums. For simplicity,
 * each thread sums an equal number of values.
 *
 * @params:
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 *      - int nthreads: The number of threads to use
 * @return:
 *      - double sum: The sum of values in 'vals' 
 * --------------------- 
 */
static double sumValues(double* vals, size_t count, int nthreads)
{
	struct thread_args thread_arg[nthreads];
	double sum = 0;

	// Print error and exit the number of values can not be split equally amongst the threads
	if (count % nthreads) {
		fprintf(stderr, "count %zu must be divisible by %d\n", count, nthreads);
		return 0;
	}

	// Split up the values in 'vals' equally and iteratively spawn the threads
	for (int i = 0; i < nthreads; i++) {
		
		// Create thread running threadSum
		//std::thread new_thread(std::ref(threadSum));

		// Initalise thread_arg structure
		thread_arg[i].count = count / nthreads;
		thread_arg[i].vals = &vals[i * count / nthreads];
		thread_arg[i].thread_sum = 0;
		//thread_arg[i].thread_ID = new_thread.get_id();
		pthread_create(&thread_arg[i].thread, NULL, threadSum, &thread_arg[i]);
	}

	// Wait for each thread to finish, then add in its partial sum to the total sum
	for (int i = 0; i < nthreads; i++) {
		//thread_arg[i].thread_ID.join();
		pthread_join(thread_arg[i].thread, NULL);
		sum += thread_arg[i].thread_sum;
	}

	return sum;
}


/*
 * Function:    allocateValues
 * ------------------------------
 * Dynamically allocates 'count' random values to an array.
 *
 * @params:
 *      - size_t count: The number of values to be created
 *      - int nthreads: The number of threads to use
 * 		- bool using_threading: If threading is used
 *  * @return:
 *      - double* vals: A pointer to the allocated values
 * --------------------- 
 */
static double* allocateValues(double count)
{
	double* vals;
	
	vals = (double*) malloc(sizeof(double) * count); // Dynamically allocate space
	if (!vals) // Exit if no space could be dynamically allocated
		exit(1);
	for (size_t i = 0; i < count; i++) { // Assign a random value to each 
		vals[i] = rand() % MAX_VALUE;
	}

	return vals;
}


/*
 * Function:    displayThreadedSumResults
 * ------------------------------
 * 
 * // Save a copy of the strings 
 *
 * @params:
 * 		- char** strings_to_write: Pointer to a list of strings
 * 								   to be written to the .csv.
 * --------------------- 
 */
static void* writeToCSV(char** strings_to_write)
{
	size_t numLines;
	std::ofstream output_file;
	
	 // Assign variables
    output_file.open(OUTPUT_FILEPATH);
	numLines = sizeof(strings_to_write); // Calculate the number of strings to write
	
	// Write strings in strings_to_write to .csv file
	for(int i = 0; i <= numLines; i++)
	{
		output_file << *strings_to_write[i]; // Write string
		free(strings_to_write[i]); // Free string inside once written
	}

	// Free memory
	free(strings_to_write);
	output_file.close();

	return NULL;


}

/*
 * Function:    getThreadedSumResults
 * ------------------------------
 * Sums the value of an array using threading. This sum is
 * then returned along with the time taken to calculate it
 * as a string, which can be printed or saved to a file.
 *
 * @params:
 * 		- double* vals: The values to be summed
 *      - size_t count: The number of values to be summed
 *      - int nthreads: The number of threads to use
 *  * @returns:
 * 		- char* results_string: String containing the results
 * --------------------- 
 */
static char* getThreadedSumResults(double* vals, size_t count, int nthreads)
{
	struct timespec start, end;
	uint64_t nanoseconds;
	double sum;
	char* results_string = (char*) malloc((MAX_STRING_SIZE + 1) * sizeof(char));

	// Calculate the sum of 'vals' and the time taken to perform this sum using threads
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
	sum = sumValues(vals, count, nthreads);
	clock_gettime(CLOCK_MONOTONIC, &end); // End timer
	nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec);

	// Print the results
	printf("thread sum: %f", sum);
	snprintf(results_string, MAX_STRING_SIZE + 1, "thread sum: %f (%d threads)\nTook %" PRIu64 
		     " ms for %zu iterations\n", sum, nthreads, nanoseconds / MS_TO_NS, count); // Write sum to a string

	return results_string;
}

/*
 * Function:    getNonThreadedSumResults
 * ------------------------------
 * Sums the value of an array without using threading. This 
 * sum is then returned along with the time taken to calculate
 * it as a string, which can be printed or saved to a file.
 *
 * @params:
 * 		- double* vals: The values to be summed
 *      - size_t count: The number of values to be summed
 * @returns:
 * 		- char* results_string: String containing the results
 * --------------------- 
 */
static char* getNonThreadedSumResults(double* vals, size_t count)
{
	struct timespec start, end;
	uint64_t nanoseconds;
	double sum;
	char* results_string = (char*) malloc((MAX_STRING_SIZE + 1) * sizeof(char));

	// Calculate the sum of 'vals' and the time taken to perform this sum without using threads
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
	sum = noThreadSum(vals, count);
	clock_gettime(CLOCK_MONOTONIC, &end); // End timer
	nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec);

	// Save results as a string
	snprintf(results_string, MAX_STRING_SIZE + 1, "no thread sum: %f\nTook %" PRIu64
	 	     " ms for %zu iterations\n", sum, nanoseconds / MS_TO_NS, count);

	return results_string;
}


/*
 * Function:    performSumming
 * ------------------------------
 * Sums the value of an array.
 * This sum is then displayed along with the time taken to
 * calculate it.
 *
 * @params:
 *      - size_t count: The number of values to be summed
 *      - int nthreads: The number of threads to use
 * 		- bool using_threading: If threading is used
 * --------------------- 
 */
static void* performSumming(size_t count, int nthreads)
{
	struct timespec start, end;
	uint64_t thread_nanoseconds;
	uint64_t no_thread_nanoseconds;
	double thread_sum;
	double no_thread_sum;
	double* vals;
	char** results = (char**) malloc((MAX_STRINGS_TO_SAVE + 1) * sizeof(char));

	// Allocate random values to be summed
	vals = allocateValues(count);

	// Calculate and save the sum of 'vals' with and without using threads
	results[THREADED_RESULT_ORDER] = getThreadedSumResults(vals, count, nthreads);
	results[NON_THREADED_RESULT_ORDER] = getNonThreadedSumResults(vals, count);
	writeToCSV(results);
	free(vals); 

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

	// Re-initalise variables based on the program's arguments
	if (argc > COUNT_ARG_POS) {
		count = atoll(argv[COUNT_ARG_POS]);
	}
	if (argc > THREADS_ARG_POS) {
		nthreads = atoi(argv[THREADS_ARG_POS]);
	}

	// Calculate the sum and print the results
	performSumming(count, nthreads); // Determine the sum using a threaded summing routine

	return 0;
}