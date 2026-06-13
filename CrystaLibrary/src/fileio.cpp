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
#include <fileio.h>

OFFFile::OFFFile(std::string file_path) {
	std::ifstream inputfile(file_path);
	std::string curr_line, token;
	std::vector<std::string> tokens;
	int points_count = 0, faces_count = 0, curr_points_count = 0, curr_faces_count=0;
	FileSections nextsection = FileSections::COUNT;
	while (!inputfile.eof()) {
		std::getline(inputfile, curr_line);
		std::istringstream tokenStream(curr_line);
		if (curr_line.empty()) continue;
		
		switch (nextsection) {
			case FileSections::COUNT: {
				// Step 1: Get number of vertices and faces
				//tokens = this->splitByDelimiter(curr_line, ' ');
				tokens = this->splitByDelimiters(curr_line, std::regex("[^\\s\\t]+"));
				points_count = std::stoi(tokens[0]);
				faces_count = std::stoi(tokens[1]);
				nextsection = FileSections::POINTLIST;
				break;
			}
			case FileSections::POINTLIST: {
				//Step 2: Get point list
				//tokens = this->splitByDelimiter(curr_line, ',');
				tokens = this->splitByDelimiters(curr_line, std::regex("[^\\s\\t,]+"));

				std::vector<double> point;
				for (std::string coordinate : tokens) {
					point.push_back(std::stod(coordinate));
				}
				this->points.push_back(point);
				++curr_points_count;
				if (curr_points_count == points_count)
					nextsection = FileSections::FACELIST;
				break;
			}
			case FileSections::FACELIST: {
				//Step 3: Get face list
				std::vector<int> face;
				//tokens = this->splitByDelimiter(curr_line, ' ');
				tokens = this->splitByDelimiters(curr_line, std::regex("[^\\s\\t]+"));

				this->indeces_counts.push_back(std::stoi(tokens[0]));
				//tokens = this->splitByDelimiter(tokens[1], ',');
				//handle the case when spaces are inbetween
				for (int i = 1; i < tokens.size(); ++i) {
					std::vector<std::string> curr_tokens = this->splitByDelimiters(tokens[i], std::regex("[^\\s\\t,]+"));
					for (std::string p_index : curr_tokens) {
						face.push_back(std::stoi(p_index));
					}
				}

				this->faces.push_back(face);
				++curr_faces_count;
				break;
			}
			default: break;
		}
	}
}

std::vector<std::string> OFFFile::splitByDelimiter(std::string s, char delimiter) {
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter)) { tokens.push_back(token); }

	return tokens;
}

std::vector<std::string> OFFFile::splitByDelimiters(std::string s, std::regex regex_delimiters) {
	std::vector<std::string> tokens;

	auto words_begin = std::sregex_iterator(s.begin(), s.end(), regex_delimiters);
	auto words_end = std::sregex_iterator();

	for (std::sregex_iterator i = words_begin; i != words_end; ++i)
		tokens.push_back((*i).str());
	return tokens;
}
