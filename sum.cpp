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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <thread>   

#define S_TO_NS 1000000000ULL 		// Conversion factor from seconds to nanoseconds
#define MS_TO_NS 1000000ULL			// Conversion factor from millseconds to nanoseconds
#define BASE_VALUES_TO_SUM 10000 	// The number of values to sum if no inputs are given to main
#define BASE_NUM_OF_THREADS 4   	// The number of threads to use if no inputs are given to main
#define COUNT_ARG_POSITION 1		// The position of the count argument ot main
#define THREADS_ARG_POSITION 2		// T

/*
 * Structure:    thread_args
 * ------------------------------
 * Defines stucture to contain the information needed
 * for arguments to the sum.cpp thread functions
 * 
 * @members:
 *      - double* vals: A pointer to the values to be summed
 *      - size_t count: The number of values to be summed
 *      - double thread_sum: The sum of values
 *      - std::thread thread: The thread used to sum values
 * ---------------------
 */
typedef struct thread_args {
	double* vals;			  
	size_t count;			   
	double thread_sum;
	std::thread thread;
} thread_args_t;



/*
 * Function:    no_thread_sum
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
static double no_thread_sum(double* vals, size_t count)
{
	double sum = 0;
	for (size_t i = 0; i < count; i++) {
		sum += vals[i];
	}
	return sum;
}

// return the sum of a single slice of values
/*
 * Function:    no_thread_sum
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
static void *thread_sum_func(void *args)
{
	struct thread_args *ta = (struct thread_args*) args;
	double sum = 0;

	for (size_t i = 0; i < ta->count; i++) {
		sum += ta->vals[i];
	}
	ta->ret = sum;
	return NULL;
}

// Break up 'vals' into a series of slices, and spawn a thread to sum each
// slice, then return the sum of these intermediate sums.
/*
 * Function:    getControlSignal
 * ------------------------------
 * Function reverses error signal for yaw and processes
 * the error signal so it will work with the method used
 * to log yaw which is from 0 to 179 and -180 to 0.
 *
 * The control signal is calculated, using PID gains and error signal
 *
 * Duty cycle limits are set for altitude and yaw so as
 * to not overload the helicopter rig and emulator.
 *
 * @params:
 *      - controller_t* piController: Pointer to the relevant
 *      conroller struct.
 *      - int32_t reference: Target yaw/altitude
 *      - int32_t measurement: Actual yaw/altitude
 *      - bool isYaw: True if the controller is for yaw. False if
 *      controller is for altitude.
 * @return:
 *      - int32_t dutyCycle: Appropriate duty cycle for the relevant
 *      PWM output as calculated by the control system
 * ---------------------
 */
static double thread_sum(double *vals, size_t count, int nthreads)
{
	// How many threads should we create?
	struct thread_args ta[nthreads];
	double sum = 0;

	// For simplicity, we're requiring them to all be equal size
	if (count % nthreads) {
		fprintf(stderr, "count %zu must be divisible by %d\n", count, nthreads);
		return 0;
	}

	// Split up the incoming data, and spawn the threads
	for (int i = 0; i < nthreads; i++) {
		ta[i].count = count / nthreads;
		ta[i].vals = &vals[i * count / nthreads];
		ta[i].ret = 0;
		if (pthread_create(&ta[i].thread, NULL, thread_sum_func, &ta[i]) < 0)
			fprintf(stderr, "Could not create thread %d\n", i);
	}

	// Wait for each thread to finish, and add in its partial sum
	for (int i = 0; i < nthreads; i++) {
		//pthread_join(ta[i].thread, NULL);
		ta[i].join();
		sum += ta[i].ret;
	}
	return sum;
}

/*
 * Function:    getControlSignal
 * ------------------------------
 * Function reverses error signal for yaw and processes
 * the error signal so it will work with the method used
 * to log yaw which is from 0 to 179 and -180 to 0.
 *
 * The control signal is calculated, using PID gains and error signal
 *
 * Duty cycle limits are set for altitude and yaw so as
 * to not overload the helicopter rig and emulator.
 *
 * @params:
 *      - controller_t* piController: Pointer to the relevant
 *      conroller struct.
 *      - int32_t reference: Target yaw/altitude
 *      - int32_t measurement: Actual yaw/altitude
 *      - bool isYaw: True if the controller is for yaw. False if
 *      controller is for altitude.
 * @return:
 *      - int32_t dutyCycle: Appropriate duty cycle for the relevant
 *      PWM output as calculated by the control system
 * ---------------------
 */
int main(int argc, char **argv)
{
	// Initalise variables
	struct timespec start, end;
	uint64_t nanoseconds;
	double *vals;
	double sum;
	size_t count = BASE_VALUES_TO_SUM;
	int nthreads = BASE_NUM_OF_THREADS;

	// Reinitalise variables based on the program's arguments
	if (argc > 1) {
		count = atoll(argv[1]);
	}
	if (argc > 2) {
		nthreads = atoi(argv[2]);
	}

	// Allocate some space for the values and initialise them with some random data
	vals = (double*) malloc(sizeof(double) * count);
	if (!vals) // Exit if no space could be dynamically allocated
		exit(1);
	for (size_t i = 0; i < count; i++) { // Assign a random value to each 
		vals[i] = rand() % 100;
	}

	// Determine the sum using a threaded summing routine
	clock_gettime(CLOCK_MONOTONIC, &start);
	sum = thread_sum(vals, count, nthreads);
	clock_gettime(CLOCK_MONOTONIC, &end);
	nanoseconds = (end.tv_sec - start.tv_sec) * SECONDS_TO_NANOSECONDS +
		(end.tv_nsec - start.tv_nsec);
	printf("thread sum: %f (%d threads)\n", sum, nthreads);
	printf("Took %" PRIu64 " ms for %zu iterations\n",
		nanoseconds / MS_TO_NS, count);

	// Determine the sum using a non-threaded summing routine
	clock_gettime(CLOCK_MONOTONIC, &start);
	sum = no_thread_sum(vals, count);
	clock_gettime(CLOCK_MONOTONIC, &end);
	nanoseconds = (end.tv_sec - start.tv_sec) * SECONDS_TO_NANOSECONDS +
		(end.tv_nsec - start.tv_nsec);
	printf("no thread sum: %f\n", sum);
	printf("Took %" PRIu64 " ms for %zu iterations\n",
		nanoseconds / MS_TO_NS, count);

	free(vals); // Free the dynamically allocated memory used to store values

	return 0;
}