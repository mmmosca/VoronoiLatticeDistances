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
#include <linearcomplex.h>

void transform_dart_to_their_dual(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, Linear_Complex_Combinatorial_Map& lc_dual_combinatorial_map, std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart)
{
	// Pick the darts handles from the Linear Complex Combinatorial map of DT and from the VORONOI
	typename Linear_Complex_Combinatorial_Map::Dart_range::iterator it1 = lc_combinatorial_map.darts().begin();
	typename Linear_Complex_Combinatorial_Map::Dart_range::iterator it2 = lc_dual_combinatorial_map.darts().begin();

	std::map<typename Linear_Complex_Combinatorial_Map::Dart_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle> dual;

	// Put all dart_handles of VORONOI in dual (DT -> VORONOI dart)
	// Dual map will have: dart_handle from DT -> dart_handle from VORONOI
	// Relate a dart in DT -> dart in VORONOI (one to one)
	for (; it1 != lc_combinatorial_map.darts().end(); ++it1, ++it2)
	{
		dual[it1] = it2;
	}

	for (typename std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>::iterator it = cell_to_dart.begin(), itend = cell_to_dart.end(); it != itend; ++it)
	{
		// first: cell handle or volume, it->second: dart_handle of LC Map
		// assoc[first]: Voronoi handle, dual[second]: dart_handle of LC Map
		// DT cell handle -> Voronoi dart_handles = LC Comb Map handle -> LC Comb Map handle 
		cell_to_dart[it->first] = dual[it->second];
	}
}

void set_geometry_of_dual(Linear_Complex_Combinatorial_Map& lc_dual_combinatorial_map, Delaunay_Triangulation& lc_combinatorial_map, std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart)
{
	for (typename std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>::iterator it = cell_to_dart.begin(), itend = cell_to_dart.end(); it != itend; ++it)
	{
		if (!lc_combinatorial_map.is_infinite(it->first)) {
			lc_dual_combinatorial_map.set_vertex_attribute(it->second, lc_dual_combinatorial_map.create_vertex_attribute(lc_combinatorial_map.dual(it->first)));
			//Point_CM p = lc_dual_combinatorial_map.point(it->second);
			//std::cout << "\nnot infinite: " << p << std::endl;
		}
		else {
			lc_dual_combinatorial_map.set_vertex_attribute(it->second, lc_dual_combinatorial_map.create_vertex_attribute());
			//Point_CM p = lc_dual_combinatorial_map.point(it->second);
			//std::cout << "\ninfinite: " << p << std::endl;
		}		
	}
}

void removeInfinteCellsFromDual(Linear_Complex_Combinatorial_Map &dual_cm, Dart_handle dh_inf) {

	// Darts belonging to Infinite volume cells 
	std::stack<Dart_handle> toremove;
	Linear_Complex_Combinatorial_Map::size_type mark_toremove = dual_cm.get_new_mark();

	// dh_inf is the starting handle in an infinite volume cell
	toremove.push(dh_inf);
	// mark the starting infinite cell
	CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 3>(dual_cm, dh_inf, mark_toremove);
	// Now we get all the volumes adjacent to the infinite volume.
	for (Linear_Complex_Combinatorial_Map::Dart_of_cell_range<3>::iterator it = dual_cm.darts_of_cell<3>(dh_inf).begin(),
		itend = dual_cm.darts_of_cell<3>(dh_inf).end(); it != itend; ++it)
	{		
		// If the next infinite cell is not marked -> mark and add it for removal
		if (!dual_cm.is_marked(dual_cm.beta(it, 3), mark_toremove)) {
			CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 3>(dual_cm, dual_cm.beta(it, 3), mark_toremove);
			toremove.push(dual_cm.beta(it, 3));
		}
	}
	// Remove all infinite cell marked
	while (!toremove.empty()) {
		dual_cm.remove_cell<3>(toremove.top());
		toremove.pop();
	}
	
	// Remove missing infinite cells
	for (Linear_Complex_Combinatorial_Map::Vertex_attribute_range::iterator v = dual_cm.vertex_attributes().begin(),
		v_end = dual_cm.vertex_attributes().end(); v != v_end; ++v)
	{
		Point_CM p = v->point();
		if ((p[0] >= INF_COMPONENT) || (p[0] <= -INF_COMPONENT) || (p[1] >= INF_COMPONENT) || (p[1] <= -INF_COMPONENT) || (p[2] >= INF_COMPONENT) || (p[2] <= -INF_COMPONENT)) {

			//std::cout << "Infinite point by verteces: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			CGAL_assertion(dual_cm.is_removable<3>(v->dart()));
			dual_cm.remove_cell<3>(v->dart());
		}
	}
	
	CGAL_assertion(dual_cm.is_without_boundary(1) && dual_cm.is_without_boundary(2));
}

void updateToFiniteSimplicialComplex(Linear_Complex_Combinatorial_Map& lc_combinatorial_map,
	Delaunay_Triangulation& DT,
	std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart) {

	// Extract the finite simplicial complex
	ACell_iterator it;
	for (it = DT.all_cells_begin(); it != DT.all_cells_end(); ++it) {
		if (DT.is_infinite(it)) {
			lc_combinatorial_map.remove_cell<3>(cell_to_dart[it]);
		}
	}
	CGAL_assertion(lc_combinatorial_map.is_valid());
}

void removeIncorrectGeometryFaces(Linear_Complex_Combinatorial_Map& lc_combinatorial_map) {
	
	std::stack<Dart_handle> toremove;
	Linear_Complex_Combinatorial_Map::size_type mark_toremove = lc_combinatorial_map.get_new_mark();

	// Remove incorrect faces within the combinatorial map
	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<2>::iterator face_dart = lc_combinatorial_map.one_dart_per_cell<2>().begin(),
		face_dart_end = lc_combinatorial_map.one_dart_per_cell<2>().end(); face_dart != face_dart_end; ++face_dart)
	{
		Dart_handle dh_edges = face_dart;
		while (!lc_combinatorial_map.is_marked(dh_edges, mark_toremove)) {
			Point_CM p = lc_combinatorial_map.point(dh_edges);
			if ((p[0] >= INF_COMPONENT) || (p[0] <= -INF_COMPONENT) || (p[1] >= INF_COMPONENT) || (p[1] <= -INF_COMPONENT) || (p[2] >= INF_COMPONENT) || (p[2] <= -INF_COMPONENT)) {
				
				//std::cout << "Infinite point by faces: " << p[0] << " " << p[1] << " " << p[2] << std::endl; ;
				CGAL_assertion(lc_combinatorial_map.is_removable<2>(face_dart));
				toremove.push(face_dart);
			}
			if (dh_edges == face_dart) {
				CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_combinatorial_map, dh_edges, mark_toremove);
			}
			dh_edges = lc_combinatorial_map.beta(dh_edges, 1);
		}
		CGAL::unmark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_combinatorial_map, dh_edges, mark_toremove);
	}

	// Remove all infinite faces
	while (!toremove.empty()) {
		lc_combinatorial_map.remove_cell<2>(toremove.top());
		toremove.pop();
	}
	
	// Remove missing infinite cells
	for (Linear_Complex_Combinatorial_Map::Vertex_attribute_range::iterator v = lc_combinatorial_map.vertex_attributes().begin(),
		v_end = lc_combinatorial_map.vertex_attributes().end(); v != v_end; ++v)
	{
		Point_CM p = v->point();
		if ((p[0] >= INF_COMPONENT) || (p[0] <= -INF_COMPONENT) || (p[1] >= INF_COMPONENT) || (p[1] <= -INF_COMPONENT) || (p[2] >= INF_COMPONENT) || (p[2] <= -INF_COMPONENT)) {

			//std::cout << "Infinite point by verteces: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			CGAL_assertion(lc_combinatorial_map.is_removable<3>(v->dart()));
			lc_combinatorial_map.remove_cell<3>(v->dart());
		}
	}

}


void display_finite_characteristics(Delaunay_Triangulation& DT) {
	std::cout << "Delaunay Triangulation ->" << 
		" 0-cells:" << DT.number_of_vertices() <<
		" 1-cells:" << DT.number_of_finite_edges() <<
		" 2-cells:" << DT.number_of_finite_facets() <<
		" 3-cells:" << DT.number_of_finite_cells() << std::endl;
}

void display_finite_characteristics(Linear_Complex_Combinatorial_Map& lc_combinatorial_map,
	Delaunay_Triangulation& DT,
	std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart) {
	
	std::map<typename Linear_Complex_Combinatorial_Map::Dart_handle, typename Delaunay_Triangulation::Cell_handle> dart_to_cell;
	std::map<typename Linear_Complex_Combinatorial_Map::Dart_handle, typename Delaunay_Triangulation::Cell_handle>::iterator it;

	for (auto& couple : cell_to_dart) {
		dart_to_cell[couple.second] = couple.first;
	}

	int count_3cells = 0, count_2cells = 0, count_1cells = 0, count_0cells = 0, notrecognized = 0, infinite_darts= 0;
	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<3>::iterator cell_dart = lc_combinatorial_map.one_dart_per_cell<3>().begin(),
		cell_dart_end = lc_combinatorial_map.one_dart_per_cell<3>().end(); cell_dart != cell_dart_end; ++cell_dart)
	{
		
		it = dart_to_cell.find(cell_dart);
		if (it != dart_to_cell.end()) {
			if (!DT.is_infinite(dart_to_cell[cell_dart])) {
				++count_3cells;
				/*
				for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<2>::iterator face_dart = lc_combinatorial_map.one_dart_per_cell<2>().begin(),
					face_dart_end = lc_combinatorial_map.one_dart_per_cell<2>().end(); face_dart != face_dart_end; ++face_dart)
				{
				}
				*/
			}
			else {
				++infinite_darts;
			}
		}
		else {
			++notrecognized;
		}		
	}

	std::cout << "Combinatorial map -> 3-cells:" << count_3cells << ", Known and Infinite darts:" << infinite_darts << ", Unknown:" << notrecognized << std::endl;
}