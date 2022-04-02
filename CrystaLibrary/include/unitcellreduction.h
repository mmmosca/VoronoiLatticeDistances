#ifndef _UNITCELLREDUCTION_H
#define _UNITCELLREDUCTION_H

//#define DEBUG
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <geom.h>

enum class UnitCellReductionMode {Niggli, Buerger};

Eigen::Matrix3d reduceUnitCell(std::vector<double> &cell_parameters, Eigen::Matrix3d &transform, bool &total_reduced);

#endif // !_UNITCELLREDUCTION_H