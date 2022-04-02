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
