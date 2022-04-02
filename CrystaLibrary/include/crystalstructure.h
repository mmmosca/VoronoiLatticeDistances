#ifndef _STRUCTURE_H
#define _STRUCTURE_H

//#define DEBUG

#include <stdio.h>
#include <Eigen/Dense>
#include <set>
#include <map>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <thread> // Standard library threads

#include <geom.h>
#include <unitcellreduction.h>
#include <geometrystructure.h>

class myUnitCell {
public:
	std::vector<double> cell_parameters;
	Eigen::Vector3d v_a, v_b, v_c;
	Eigen::Matrix3d fractToCart;
		
	double getVolume();
};


class Lattice : public GeometryStructure, public myUnitCell {
private:	
	//	Update the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	void updateCellParameters();
		
	//	Update and reduce the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	bool updateReducedCellParameters();

public:
	Lattice() : GeometryStructure() {};

	std::vector<double> getCellParameters();

	void setCellParameters(std::vector<double> cell_params);

	//	Update the Lattice to Niggli reduced mode
	bool updateParametersToReducedMode();

	/*
	*	Update the Lattice to Original mode, only cell vectors are updated.
	*	Use spanTheLattice to recalculate points
	*/
	void updateParametersToOrginalMode();

	/*
	*	Generate points in the Lattice around the origin by n copies in all three dimensions
	*	param n: number of unit cells for each direction
	*	param positive: If true unit cell will be extended towards positive directions
	*/
	void spanTheLattice(int n, bool positive = false);

	// Clear the instance of the Lattice
	void clearTheLattice();

	//Print info
	void print_lattice_info();
};

#endif // !_STRUCTURE_H