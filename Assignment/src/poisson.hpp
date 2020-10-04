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
    private:
        int allocateUserInputs(int, char**);
        void* allocateVolume(void);
} poisson_args_t;


/*
 * Function:    poisson_dirichlet
 * ------------------------------
 * Solve Poisson's equation for a rectangular box with Dirichlet
 * boundary conditions on each face. The remaining documentation
 * is included to fulfill the ENCE464 Assignment 2 specifications.
 * However, future program develpments would have this function
 * as a method to the poisson_args_t object.
 * 
 * --------------------- 
*/
/// Solve Poisson's equation for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param Vbound is the potential on the boundary
/// \param xsize is the number of elements in the x-direction
/// \param ysize is the number of elements in the y-direction
/// \param zsize is the number of elements in the z-direction
/// \param delta is the voxel spacing in all directions
/// \param numiters is the number of iterations to perform
/// \param numcores is the number of CPU cores to use.  If 0, an optimal number is chosen
void poisson_dirichlet (double * __restrict__ source,
                        double * __restrict__ potential,
                        double Vbound,
                        unsigned int xsize, unsigned int ysize, unsigned int zsize, double delta,
                        unsigned int numiters, unsigned int numcores);
                        

#endif /* POISSON_H */
