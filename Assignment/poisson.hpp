/* ****************************************************************
 * poisson.hpp
 *
 * Part of ENCE464 Assignment 2
 * Header file containing the classes and functions needed to
 * calculate the electric potential 3-D cube, for a given charge
 * distribution. This is solved through Poisson's equation using
 * Jacobi relaxation. 
 *
 * Based off poisson.hpp - Michael Hayes, UC ECE
 * 
 * ENCE464 Assignment 2 Group 1
 * Creators: Matt Blake          58979250
 *           Derrick Edward      18017758
 * Last modified: 03/10/2020
 * ***************************************************************/

#ifndef POISSON_H
#define POISSON_H


/*
 * Class:    poisson_args
 * ------------------------------
 * Contains the information needed to solve Poisson's
 * equation.
 * 
 * @members:
        - unsigned int num_iters: Number of iterations of Jacobi relaxation to perform.
        - unsigned int num_cores: Number of CPU cores to use. If 0/undefined, an
                                  optimal number is chosen.
        - unsigned int x_size: Size of the rectangular box used in the x dimension.
        - unsigned int y_size: Size of the rectangular box used in the y dimension.
        - unsigned int z_size: Size of the rectangular box used in the z dimension .   
        - double delta: The spacing (in all directions) between voxels (meters).
        - double V_bound: The boundary condition potential.
    @methods:
        - initPoissonArgs() - Initialises a poisson_args object
        - allocateUserInputs() - Allocates the user defined members during initalisation.
        - allocateVolume() - Allocates memory for the 3-D cuboid's volume.
        - poissonDirichlet() - Solves Poisson's equation using Jacobi relaxation.
        - performJacobiIteration() - Performs an iteration of Jacobi relaxation.
 * ---------------------
 */
typedef class poisson_args
{
	public:
        unsigned int num_iters;
        unsigned int num_cores;    
        unsigned int x_size;
        unsigned int y_size;
        unsigned int z_size;    
        double delta;
        double V_bound;
        double* source;
        double* potential;
        void* initPoissonArgs(int, char**);
        int poissonDirichlet(void);
    private:
        int allocateUserInputs(int, char**);
        void* allocateVolume(void);
        void* performJacobiIteration(double*);
} poisson_args_t;


#endif /* POISSON_H */
