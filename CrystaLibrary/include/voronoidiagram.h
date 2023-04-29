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
