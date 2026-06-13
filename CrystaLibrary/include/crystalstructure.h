/*Copyright (C) 2018-2022,2026 Marco M. Mosca

This file is part of VoronoiLatticeDistances.

VoronoiLatticeDistances is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VoronoiLatticeDistances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VoronoiLatticeDistances. If not, see <https://www.gnu.org/licenses/>.
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
};

#endif // !_STRUCTURE_H
