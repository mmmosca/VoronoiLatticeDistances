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
#ifndef _EIGENSUPPORT_H
#define _EIGENSUPPORT_H

struct EigenVCompare {
	bool operator() (const Eigen::VectorXd& u, const Eigen::VectorXd& v) const {
		bool cond = false;
		for (int i = 0; i < u.size(); ++i) {
			if (u(i) < v(i)) {
				cond = true;
				break;
			}
			else if (u(i) > v(i)) {
				cond = false;
				break;
			}
		}
		return cond;
	}
};

#endif // !_EIGENSUPPORT_H
