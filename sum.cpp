#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include <thread>   

// structure we're going to use for arguments to our thread functions
struct thread_args {
	double *vals;
	size_t count;
	double ret;
	std::thread thread;
};

// no_thread_sum returns the sum of the values in 'vals', of length 'count'
static double no_thread_sum(double *vals, size_t count)
{
	double sum = 0;
	for (size_t i = 0; i < count; i++) {
		sum += vals[i];
	}
	return sum;
}

// return the sum of a single slice of values
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


int main(int argc, char **argv)
{
	struct timespec start, end;
	size_t count = 10000;
	double *vals;
	double sum;
	uint64_t nanoseconds;
	int nthreads = 4;

	if (argc > 1) {
		count = atoll(argv[1]);
	}
	if (argc > 2) {
		nthreads = atoi(argv[2]);
	}

	/* Allocate some space for our values, and initialise them with some
	 * random data
	 */
	vals = (double*) malloc(sizeof(double) * count);
	if (!vals)
		exit(1);

	for (size_t i = 0; i < count; i++) {
		vals[i] = rand() % 100;
	}

	// Determine the sum using a threaded summing routine
	clock_gettime(CLOCK_MONOTONIC, &start);
	sum = thread_sum(vals, count, nthreads);
	clock_gettime(CLOCK_MONOTONIC, &end);
	nanoseconds = (end.tv_sec - start.tv_sec) * 1000000000ULL +
		(end.tv_nsec - start.tv_nsec);
	printf("thread sum: %f (%d threads)\n", sum, nthreads);
	printf("Took %" PRIu64 " ms for %zu iterations\n",
		nanoseconds / 1000000, count);

	// Determine the sum using a non-threaded summing routine
	clock_gettime(CLOCK_MONOTONIC, &start);
	sum = no_thread_sum(vals, count);
	clock_gettime(CLOCK_MONOTONIC, &end);
	nanoseconds = (end.tv_sec - start.tv_sec) * 1000000000ULL +
		(end.tv_nsec - start.tv_nsec);
	printf("no thread sum: %f\n", sum);
	printf("Took %" PRIu64 " ms for %zu iterations\n",
		nanoseconds / 1000000, count);

	free(vals);

	return 0;
}