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
#include <cifio.h>

/**************************\
|*******  PUBLIC  *********|
\**************************/
const std::string CIFHandler::cifexample = "cifexample.cif";

const std::vector<std::string> CIFHandler::cell_parameters_tags = { "_cell_length_a", "_cell_length_b", "_cell_length_c",
							"_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma" };

int CIFHandler::findValidBlockIndex(gemmi::cif::Document doc) {
	for (int i = 0; i < doc.blocks.size(); ++i) {
		// Check if it is the right block
		if (doc.blocks[i].find("_atom_site_", { "label", "fract_x", "fract_y" , "fract_z" }).get_loop()) {
			return i;
		}
	}
	return -1;
}
