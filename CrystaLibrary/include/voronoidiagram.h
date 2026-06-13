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
#ifndef _VORONOIDIAGRAM_H
#define _VORONOIDIAGRAM_H

#include <voronoicell.h>

class VoronoiDiagram {
public:
	std::vector<VoronoiCell> voronoicells_vector;

	VoronoiDiagram(Lattice& lattice, int extend);
	std::vector<VoronoiCell> getVoronoiCellsFromPointsInside(std::vector<Eigen::VectorXd> points_inside);
	std::vector<VoronoiCell> getVoronoiCellsFromPointsInside_Ord(std::vector<Eigen::VectorXd> points_inside);
	void writeVoronoiDomainsToVTP(std::string filename);
	void writeVoronoiDomainsToOFF(std::string filename);

};
#endif // !_VORONOIDIAGRAM_H
