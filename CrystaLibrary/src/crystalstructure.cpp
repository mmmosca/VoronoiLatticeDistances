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
#include "crystalstructure.h"

/**************************\
|*******  UNITCELL  *******|
\**************************/

double myUnitCell::getVolume() {
	return abs((Eigen::Vector3d(this->v_a).cross(Eigen::Vector3d(this->v_b))).dot(this->v_c));
}


/**************************\
|*******  LATTICE  ********|
\**************************/

std::vector<double> Lattice::getCellParameters() {
	return this->cell_parameters;
}

void Lattice::setCellParameters(std::vector<double> cell_params) {
	this->cell_parameters.assign(cell_params.begin(), cell_params.end());
}

void Lattice::updateCellParameters() {
	Eigen::Matrix3d axis;

	assert(this->cell_parameters.size() == 6);
	axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	this->v_a = axis.col(0);
	this->v_b = axis.col(1);
	this->v_c = axis.col(2);
}

bool Lattice::updateReducedCellParameters() {
	Eigen::Matrix3d fixedX_axis, reduced_transform, reduced_axis;
	bool total_reduced = false;
	
	assert(this->cell_parameters.size() == 6);
	reduced_axis = reduceUnitCell(this->cell_parameters, reduced_transform, total_reduced);
	// The following conversion is needed to fix A as X axis
	fixedX_axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	this->v_a = fixedX_axis.col(0);
	this->v_b = fixedX_axis.col(1);
	this->v_c = fixedX_axis.col(2);

	return total_reduced;
}

bool Lattice::updateParametersToReducedMode() {
	this->clearCombinations();
	return this->updateReducedCellParameters();

}

void Lattice::updateParametersToOrginalMode() {
	this->clearCombinations();
	this->updateCellParameters();
}

void Lattice::spanTheLattice(int n, bool positive) {
	//assert(n > 1);
	this->clearCombinations();
	this->addGenerator(this->v_a);
	this->addGenerator(this->v_b);
	this->addGenerator(this->v_c);
	if (positive) {
		this->updateCoefficientsForPositiveDirection(n);
	}
	else {
		// n = i * 2 + 1
		// i = (n-1)/2 = floor(n/2)
		// For n = 2, I want coefficients [-1,0,1]
		// n = 3, I want coefficients [-1,0,1,2]
		// n = 4, I want coefficients [-2,-1,0,1,2]
		// n = 5, [-2,-1,0,1,2,3]
		int k = n / 2;
		if (n % 2 == 0) {
			this->updateCoefficientsForEveryDirection(k);
		}
		else {
			this->updateCoefficientsFromRange(-k,k+1);
		}
	}
	this->updateCoefficientCombinations();
	this->updateGeneratorCombinations();
}

void Lattice::clearTheLattice() {
	this->clearCombinations();
	this->cell_parameters.clear();
	this->v_a.setZero();
	this->v_b.setZero();
	this->v_c.setZero();
}
