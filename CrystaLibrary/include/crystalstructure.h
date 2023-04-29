/*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*/
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