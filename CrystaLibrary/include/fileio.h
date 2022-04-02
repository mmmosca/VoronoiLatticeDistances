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
