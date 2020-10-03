/* ****************************************************************
 * cache_example.cpp
 *
 * Shows linear memory access, and non-linear memory access.
 * Designed to be used with cachegrind to show read cache misses.
 *
 * Based off cache_grind.c - Andre Renaud, UC ECE Department, 2020
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 27/09/2020
 *
 * ***************************************************************/

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>



/**
 * Sum a matrix of numbers with linear memory access
 */
int iterate_good(int *buffer, size_t x_max, size_t y_max)
{
	int sum = 0;
	printf("Running good iteration\n");
	for (size_t y = 0; y < y_max; y++) {
		for (size_t x = 0; x < x_max; x++) {
			sum = buffer[y * x_max + x];
		}
	}
	return sum;
}



/**
 * Sum a matrix of numbers with non-linear memory access
 */
int iterate_bad(int *buffer, size_t x_max, size_t y_max)
{
	int sum = 0;
	printf("Running bad iteration\n");
	for (size_t x = 0; x < x_max; x++) {
		for (size_t y = 0; y < y_max; y++) {
			sum = buffer[y * x_max + x];
		}
	}
	return sum;
}



/*
 * Function:    main
 * ------------------------------
 * The main function of cache_example.cpp. This function sums a 
 * buffer of ints using either linear or non-linear memory access.b
 *
 * @params:
 *      - int argc: The number of arguments passed into main
 *      - char** argv: The vector of arguments passed into main.
 * 					   This should contain if 'good' (linear) or
 * 					   'bad' (non-linear) memory access should be
 * 						used.
 * ---------------------
 */

int main(int argc, char **argv)
{
	size_t x_max = 5000;
	size_t y_max = 5000;
	int* buffer = (int *) calloc(x_max * y_max, sizeof(int));
	int good = 0;
	int sum;

	if (argc >= 2) {
		good = strcmp(argv[1], "good") == 0;
	}
	if (good) { // Calculate sum of buffer using linear methods
		sum = iterate_good(buffer, x_max, y_max);
	} else {  // Calculate sum of buffer using non-linear methods
		sum = iterate_bad(buffer, x_max, y_max);
	}

	std::cout << sum << '\n'; // Print the sum to the terminal

	return 0;
}
