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
#ifndef _FILEIO_H
#define _FILEIO_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>
#include <iostream>

class OFFFile {
private:
	enum class FileSections {COUNT, POINTLIST, FACELIST};
	std::vector<std::vector<double>> points;
	std::vector<std::vector<int>> faces;
	std::vector<int> indeces_counts;

	std::vector<std::string> splitByDelimiter(std::string s, char delimiter);
	std::vector<std::string> splitByDelimiters(std::string s, std::regex regex_delimiters);
public:
	OFFFile(std::string file_path);
	std::vector<std::vector<double>> getPoints() { return this->points; }
	std::vector<std::vector<int>> getFaces() { return this->faces; }
};

#endif // !_FILEIO_H
