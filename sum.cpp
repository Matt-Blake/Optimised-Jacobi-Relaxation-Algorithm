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
#include <thread>   

#define S_TO_NS 			1000000000ULL 	// Conversion factor from seconds to nanoseconds
#define MS_TO_NS 			1000000ULL		// Conversion factor from millseconds to nanoseconds
#define BASE_VALUES_TO_SUM 	10000 			// The number of values to sum if no inputs are given to main
#define BASE_NUM_THREADS 	4   			// The number of threads to use if no inputs are given to main
#define COUNT_ARG_POS 		1				// The position of the count argument to main
#define THREADS_ARG_POS 	2				// The position of the thread argument to main
#define MAX_VALUE			100				// The values will being summed will be between 0 and (MAX_VALUE - 1)
#define USING_THREADS		true 			// Boolean defining if threads are being used
#define NOT_USING_THREADS	false			// Boolean defining if threads are not being used
#define EXIT_STATUS			1				// The status value returned when the program is exited due to a malloc error


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
 *      - std::thread thread: The thread used to sum values
 * ---------------------
 */
typedef struct thread_args
{
	double* vals;			  
	size_t count;			   
	double thread_sum;
	std::thread thread;
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

		// Initalise thread_arg_t structure
		thread_arg[i].count = count / nthreads;
		thread_arg[i].vals = &vals[i * count / nthreads];
		thread_arg[i].thread_sum = 0;
		thread_arg[i].thread = thread();
	}

	// Wait for each thread to finish, then add in its partial sum to the total sum
	for (int i = 0; i < nthreads; i++) {
		thread_arg[i].thread.join();
		sum += thread_arg[i].thread_sum;
	}

	return sum;
}


/*
 * Function:    displaySumResults
 * ------------------------------
 * Sums the value of an array (with or without using threading).
 * This sum is then displayed along with the time taken to
 * calculate it.
 *
 * @params:
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 *      - int nthreads: The number of threads to use
 * 		- bool using_threading: If threading is used
 * --------------------- 
 */
static double displaySumResults(double* vals, size_t count, int nthreads, bool using_threading)
{
	struct timespec start, end;
	uint64_t nanoseconds;
	double sum;

	// Calculate the sum of 'vals' and the time taken to perform this sum
	clock_gettime(CLOCK_MONOTONIC, &start); // Start timer
	if (using_threading) {
		sum = sumValues(vals, count, nthreads);
	} else {
		noThreadSum(vals, count);
	}
	clock_gettime(CLOCK_MONOTONIC, &end); // End timer
	nanoseconds = (end.tv_sec - start.tv_sec) * S_TO_NS + (end.tv_nsec - start.tv_nsec);

	// Print results to the output
	if (using_threading) {
		printf("thread sum: %f (%d threads)\n", sum, nthreads);
	} else {
		printf("no thread sum: %f\n", sum);
	}
	printf("Took %" PRIu64 " ms for %zu iterations\n", nanoseconds / MS_TO_NS, count);
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
		vals[i] = rand() % 100;
	}

	return vals;
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
	double *vals;
	size_t count = BASE_VALUES_TO_SUM;
	int nthreads = BASE_NUM_THREADS;

	// Re-initalise variables based on the program's arguments
	if (argc > COUNT_ARG_POS) {
		count = atoll(argv[COUNT_ARG_POS]);
	}
	if (argc > THREADS_ARG_POS) {
		nthreads = atoi(argv[THREADS_ARG_POS]);
	}
	vals = allocateValues(count); // Allocate random values to be summed

	// Calculate the sum and print the results
	displaySumResults(vals, count, nthreads, USING_THREADS); // Determine the sum using a threaded summing routine
	displaySumResults(vals, count, nthreads, NOT_USING_THREADS); // Determine the sum without using threads

	// Free the dynamically allocated memory used to store values
	free(vals); 

	return 0;
}