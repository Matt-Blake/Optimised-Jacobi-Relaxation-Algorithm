// poisson.hpp

#ifndef POISSON_H
#define POISSON_H

// Solve Poisson's equation for a rectangular box with Dirichlet
// boundary conditions on each face.
void poissonDirichlet (double *__restrict__ source,
                        double *__restrict__ potential,
                        double V_bound,
                        unsigned int x_size, unsigned int y_size, unsigned int z_size,
                        double delta, unsigned int max_iters, unsigned int num_cores);
#endif /* POISSON_H */
