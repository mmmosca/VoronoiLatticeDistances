#ifndef _PERIODICPOINTCLOUD_H
#define _PERIODICPOINTCLOUD_H

#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>
#include <eigensupport.h>

/*
*	Periodic Point Cloud class to store basis and their linear combinations.
*	Coefficients for linear combinations are integers.
*	This class handles the Periodic structure in Cartesian or Fractional coordinates.
*/
class PeriodicPointCloud {

public:

	// Set of unique generators
	std::set<Eigen::VectorXd, EigenVCompare> G;

	// Set of generators combinations
	std::vector<Eigen::VectorXd> G_combinations;

	// Vector of integer coefficients
	std::vector<int> C;

	// Set of coefficients permutations
	std::vector<Eigen::VectorXd> C_combinations;

	// Vectors dimension, number of generators, number of linear combinations, number of coefficients permutaions, number of coefficients
	int dimension, generator_count = 0, g_combination_count = 0, c_combination_count = 0, coef_count = 0;

	// Param d: generators dimension
	PeriodicPointCloud(const int& d) : dimension{ d } {};

	std::vector<Eigen::VectorXd>& getPointCloud();

	// Generate all the combinations with the available generators and coefficients
	void updateGeneratorCombinations();

	void updateGeneratorCombinationsForPoints(std::vector<Eigen::VectorXd> points);

	void extendBasisVectorsCoefficientsAndLinearCombinations(int n, std::vector<Eigen::VectorXd> points = {});

	/*
	*	Add coefficients for combinations towards every direction
	*	n is the number of positive integer coefficients (negatives are added as well [-n, n])
	*/
	void updateCoefficientsForEveryDirection(int n);

	/*
	*	Add coefficients for combinations towards only positive directions
	*	n is the number of positive integer coefficients [0, n])
	*/
	void updateCoefficientsForPositiveDirection(int n);

	void updateCoefficientsFromRange(int i, int j);

	//	Update the coefficient combinations to be used for generator combinations
	void updateCoefficientCombinations(std::vector<int> to_filter = {}, Eigen::VectorXd v = {}, int n = -1);

	//	Add basis in the structure
	void addGenerator(Eigen::VectorXd g);

	// Clear all the structure
	void clear();

	void clearCoefficients();

	// Clear just the linear combinations
	void clearCombinations();

	//	Print some info about the PeriodicPointCloud class
	void print_info();
};

#endif // !_PERIODICPOINTCLOUD_H