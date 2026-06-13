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
#ifndef _UNITCELLREDUCTION_H
#define _UNITCELLREDUCTION_H

#include "logging.h"
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <geom.h>

enum class UnitCellReductionMode {Niggli, Buerger};

Eigen::Matrix3d reduceUnitCell(std::vector<double> &cell_parameters, Eigen::Matrix3d &transform, bool &total_reduced);

#endif // !_UNITCELLREDUCTION_H
