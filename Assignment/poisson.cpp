#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/// Solve Poisson's equation for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param V_bound is the potential on the boundary
/// \param x_size is the number of elements in the x-direction
/// \param y_size is the number of elements in the y-direction
/// \param z_size is the number of elements in the z-direction
/// \param delta is the voxel spacing in all directions
/// \param num_iters is the number of iterations to perform
/// \param num_cores is the number of CPU cores to use.  If 0, an optimal number is chosen
void poissonDirichlet (double * __restrict__ source,
                        double * __restrict__ potential,
                        double V_bound,
                        unsigned int x_size, unsigned int y_size, unsigned int z_size, double delta,
                        unsigned int num_iters, unsigned int num_cores)
{
    // source[i, j, k] is accessed with source[((k * ysize) + j) * xsize + i]
    // potential[i, j, k] is accessed with potential[((k * ysize) + j) * xsize + i]    
    size_t size = (size_t)y_size * z_size * x_size * sizeof(double);
	double *input = (double *)malloc(size);
	if (!input) {
		fprintf(stderr, "malloc failure\n");
		return;
	}
	memcpy(input, source, size);
	for (unsigned int iter = 0; iter < num_iters; iter++) {
		for (unsigned int x = 0; x < x_size; x++) {
			for (unsigned int z = 0; z < z_size; z++) {
				for (unsigned int y = 0; y < y_size; y++) {
					double res = 0;

					if (x < x_size - 1)
						res += input[((z * y_size) + y) * x_size + (x + 1)];
					else
						res += V_bound;
					if (x > 0)
						res += input[((z * y_size) + y) * x_size + (x - 1)];
					else
						res += V_bound;

					if (y < y_size - 1)
						res += input[((z * y_size) + (y + 1)) * x_size + x];
					else
						res += V_bound;
					if (y > 0)
						res += input[((z * y_size) + (y - 1)) * x_size + x];
					else
						res += V_bound;

					if (z < z_size - 1)
						res += input[(((z + 1) * y_size) + y) * x_size + x];
					else
						res += V_bound;
					if (z > 0)
						res += input[(((z - 1) * y_size) + y) * x_size + x];
					else
						res += V_bound; 

					res -= delta * delta * source[((z * y_size) + y) * x_size + x];

					res /= 6;

					potential[((z * y_size) + y) * x_size + x] = res;
				}
			}
		}
		memcpy(input, potential, size);
	}
	free(input);
}
