/*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
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
	void printFile();
};

#endif // !_FILEIO_H
