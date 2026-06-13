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
#ifndef _CIF_IO_H
#define _CIF_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>

#include "gemmi/cif.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/symmetry.hpp"

class CIFHandler {
private:

	const static std::string cifexample;

public:

	const static std::vector<std::string> cell_parameters_tags;

	static int findValidBlockIndex(gemmi::cif::Document doc);

};

#endif // !_CIF_IO_H
