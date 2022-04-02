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
