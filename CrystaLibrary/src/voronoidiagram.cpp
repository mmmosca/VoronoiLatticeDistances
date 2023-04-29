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
#include <voronoidiagram.h>

/********************************\
|******* VORONOIDIAGRAM *********|
\********************************/
VoronoiDiagram::VoronoiDiagram(Lattice& lattice, int extend)
{
	lattice.updateParametersToReducedMode();
	lattice.spanTheLattice(extend);

	std::vector<std::vector<Point_CM>> polyhedra_points;
	std::vector<std::vector<std::vector<int>>> polyhedra_faces;
	Linear_Complex_Combinatorial_Map voronoi_lc_cm = lattice.getFiniteVoronoiDiagram();
	//GeometryStructure::getPolyhedraData_New(voronoi_lc_cm, polyhedra_points, polyhedra_faces, motif.limit_inf, motif.limit_sup);
	GeometryStructure::getPolyhedraData(voronoi_lc_cm, polyhedra_points, polyhedra_faces);
	// Assert the number of volumes is valid
	int vcell_count = polyhedra_points.size();
	assert(vcell_count == polyhedra_faces.size());
	// Iterate vectors of polyhedra
	// Make VoronoiCell instances
	for (int i = 0; i < vcell_count; ++i) {
		VoronoiCell vcell;
		vcell.setPoints(polyhedra_points[i]);
		vcell.setFaces(polyhedra_faces[i]);
		vcell.fixGeometry();
		if (!vcell.hasCorrectGeometry()) continue;
		vcell.updatePolyhedron();
		this->voronoicells_vector.push_back(vcell);
		//std::cout << "\nPolyhedra count: " << vcell_count <<" -- Polyhedron points' count:" << vcell.getPoints().size() << " -- Polyhedron faces' count: " << vcell.getFaces().size() << "\n";
	}
}

std::vector<VoronoiCell> VoronoiDiagram::getVoronoiCellsFromPointsInside(std::vector<Eigen::VectorXd> points_inside) {
	std::vector<VoronoiCell> queried_voronois(points_inside.size());
	for (int i = 0; i < this->voronoicells_vector.size(); ++i) {
		if (this->voronoicells_vector[i].getPoints().empty()) throw std::runtime_error{"Empty Voronoi Domain"};
		std::vector<int> point_indeces = this->voronoicells_vector[i].getPointIndecesInside(points_inside);
		// There must be 1 point inside a voronoi domain, otherwise incorrect geometry
		if (point_indeces.size() == 1) {
			this->voronoicells_vector[i].setSite(points_inside[point_indeces[0]]);
			queried_voronois[point_indeces[0]] = this->voronoicells_vector[i];
		}
	}
	for (int i = 0; i < queried_voronois.size(); ++i) {
		if (queried_voronois[i].getPoints().empty()) {
			queried_voronois.erase(queried_voronois.begin() + i);
		}
	}
	return queried_voronois;
}

std::vector<VoronoiCell> VoronoiDiagram::getVoronoiCellsFromPointsInside_Ord(std::vector<Eigen::VectorXd> points_inside) {
	std::vector<VoronoiCell> queried_voronois;
	for (int j = 0; j < points_inside.size(); ++j) {
		for (int i = 0; i < this->voronoicells_vector.size(); ++i) {
			// check on size of points and faces already done in the constructor
			// If point is inside, voronoi found -> break
			if (this->voronoicells_vector[i].isPointInside(points_inside[j])) {
				this->voronoicells_vector[i].setSite(points_inside[j]);
				queried_voronois.push_back(this->voronoicells_vector[i]);
				break;
			}
		}
	}
	return queried_voronois;
}

void VoronoiDiagram::writeVoronoiDomainsToVTP(std::string filename)
{
	for (int i = 0; i < this->voronoicells_vector.size(); ++i) {
		std::string vcell_name(filename + "_" + std::to_string(i) + ".vtp");
		this->voronoicells_vector[i].writeToVTPFileFormat(vcell_name.c_str());
	}
}

void VoronoiDiagram::writeVoronoiDomainsToOFF(std::string filename)
{
	for (int i = 0; i < this->voronoicells_vector.size(); ++i) {
		std::string vcell_name(filename + "_" + std::to_string(i) + ".off");
		this->voronoicells_vector[i].writeToOFFFileFormat(vcell_name.c_str());
	}
}