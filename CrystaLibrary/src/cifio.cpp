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
