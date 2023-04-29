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
#include <geometrystructure.h>


/*************************************\
|*******  GEOMETRY STRUCTURE  ********|
\*************************************/

std::vector<Eigen::VectorXd> GeometryStructure::getPoints() {
	return this->G_combinations;
}

Linear_Complex_Combinatorial_Map GeometryStructure::getDelaunayTriangulation() {
	// Compute Delaunay triangulation
	Delaunay_Triangulation DT;
	// Update points from the lattice to Delaunay Points
	std::vector<Point_CM> dl_points;

	// Add point of the lattice
	for (auto &p : this->G_combinations) {
		Point_CM p_cm(p[0], p[1], p[2]);
		dl_points.push_back(p_cm);
	}

	DT.insert(dl_points.begin(), dl_points.end());
	CGAL_assertion(DT.is_valid(false));
	// Convert the triangulation to a Linear Complex Combinatorial Map
	Linear_Complex_Combinatorial_Map cm;
	// Make a map from the cell (3-simplex : volume) -> dart handle (all the volumes have a dart handle)
	std::map<Delaunay_Triangulation::Cell_handle, Linear_Complex_Combinatorial_Map::Dart_handle> cell_to_dart;
	// Update the combinatorial map: Tetrahedra are added to the Linear Complex combinatorial map (No change in darts). Returns a dart handle to an infinite cell (3-simplex)
	Dart_handle dh = CGAL::import_from_triangulation_3<Linear_Complex_Combinatorial_Map, Delaunay_Triangulation>(cm, DT, &cell_to_dart);

	//std::cout << "Delaunay triangulation of " << dl_points.size() << " vertices with infinite vertex" << std::endl;
	//cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << cm.is_without_boundary() << ", valid=" << cm.is_valid() << std::endl;;

	return cm;
}

Linear_Complex_Combinatorial_Map GeometryStructure::getFiniteDelaunayTriangulation() {
	// Compute Delaunay triangulation
	Delaunay_Triangulation DT;
	// Update points from the lattice to Delaunay Points
	std::vector<Point_CM> dl_points;

	// Add point of the lattice
	for (auto &p : this->G_combinations) {
		Point_CM p_cm(p[0], p[1], p[2]);
		dl_points.push_back(p_cm);
	}

	DT.insert(dl_points.begin(), dl_points.end());
	CGAL_assertion(DT.is_valid(false));
	// Convert the triangulation to a Linear Complex Combinatorial Map
	Linear_Complex_Combinatorial_Map cm;
	// Make a map from the cell (3-simplex : volume) -> dart handle (all the volumes have a dart handle)
	std::map<Delaunay_Triangulation::Cell_handle, Linear_Complex_Combinatorial_Map::Dart_handle> cell_to_dart;
	// Update the combinatorial map: Tetrahedra are added to the Linear Complex combinatorial map (No change in darts). Returns a dart handle to an infinite cell (3-simplex)
	Dart_handle dh = CGAL::import_from_triangulation_3<Linear_Complex_Combinatorial_Map, Delaunay_Triangulation>(cm, DT, &cell_to_dart);

	//std::cout << "Finite Delaunay triangulation of " << dl_points.size() << " vertices without infinite vertex" << std::endl;
	//cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << cm.is_without_boundary() << ", valid=" << cm.is_valid() << std::endl;
	//display_finite_characteristics(DT);
	//display_finite_characteristics(cm, DT, cell_to_dart);
	//updateToFiniteSimplicialComplex(cm, DT, cell_to_dart);
	removeIncorrectGeometryFaces(cm);

	return cm;
}

Linear_Complex_Combinatorial_Map GeometryStructure::getVoronoiDiagram(Dart_handle &dh_inf) {
	// Compute Delaunay triangulation
	Delaunay_Triangulation DT;
	// Update points from the lattice to Delaunay Points
	std::vector<Point_CM> dl_points;

	// Add point of the lattice
	for (auto &p : this->G_combinations) {
		Point_CM p_cm(p[0], p[1], p[2]);
		dl_points.push_back(p_cm);
	}

	DT.insert(dl_points.begin(), dl_points.end());

	CGAL_assertion(DT.is_valid(false));

	// Convert the triangulation to a Linear Complex Combinatorial Map
	Linear_Complex_Combinatorial_Map cm;
	// Make a map from the cell (3-simplex : volume) -> dart handle (all the volumes have a dart handle)
	std::map<Delaunay_Triangulation::Cell_handle, Linear_Complex_Combinatorial_Map::Dart_handle> cell_to_dart;
	// Update the combinatorial map: Tetrahedra are added to the Linear Complex combinatorial map (No change in darts). Returns a dart handle to an infinite cell (3-simplex)
	Dart_handle dh = CGAL::import_from_triangulation_3<Linear_Complex_Combinatorial_Map, Delaunay_Triangulation>(cm, DT, &cell_to_dart);

	//std::cout << "Delaunay triangulation of " << dl_points.size() << " vertices with " << DT.number_of_finite_cells() << " number of finite cell" << std::endl;
	//cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << cm.is_without_boundary() << ", valid=" << cm.is_valid() << std::endl;;
	//cm.display_darts(std::cout);

	// Compute the Voronoi from Delaunay triangulation
	Linear_Complex_Combinatorial_Map dual_cm;
	// Compute the dual of the Combinatorial map cm
	Dart_handle dh_dual = cm.dual(dual_cm, dh);

	//Update the darts and vertices of Voronoi cells
	transform_dart_to_their_dual(cm, dual_cm, cell_to_dart);
	set_geometry_of_dual(dual_cm, DT, cell_to_dart);

	dh_inf = dh_dual;
	//std::cout << "Voronoi Diagram of " << dl_points.size() << " vertices with infinite vertex" << std::endl;
	//dual_cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << dual_cm.is_without_boundary() << ", valid=" << dual_cm.is_valid() << std::endl;;

	return dual_cm;
}

Linear_Complex_Combinatorial_Map GeometryStructure::getFiniteVoronoiDiagram() {
	// Compute Delaunay triangulation
	Delaunay_Triangulation DT;
	// Update points from the lattice to Delaunay Points
	std::vector<Point_CM> dl_points;

	// Add point of the lattice
	for (auto &p : this->G_combinations) {
		Point_CM p_cm(p[0], p[1], p[2]);
		dl_points.push_back(p_cm);
	}

	DT.insert(dl_points.begin(), dl_points.end());
	CGAL_assertion(DT.is_valid(false));

	// Convert the triangulation to a Linear Complex Combinatorial Map
	Linear_Complex_Combinatorial_Map cm;
	// Make a map from the cell (3-simplex : volume) -> dart handle (all the volumes have a dart handle)
	std::map<Delaunay_Triangulation::Cell_handle, Linear_Complex_Combinatorial_Map::Dart_handle> cell_to_dart;
	// Update the combinatorial map: Tetrahedra are added to the Linear Complex combinatorial map (No change in darts). Returns a dart handle to an infinite cell (3-simplex)
	Dart_handle dh = CGAL::import_from_triangulation_3<Linear_Complex_Combinatorial_Map, Delaunay_Triangulation>(cm, DT, &cell_to_dart);

	//std::cout << "Delaunay triangulation of " << dl_points.size() << " vertices with " << DT.number_of_finite_cells() << " number of finite cell" << std::endl;
	//cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << cm.is_without_boundary() << ", valid=" << cm.is_valid() << std::endl;;

	// Compute the Voronoi from Delaunay triangulation
	Linear_Complex_Combinatorial_Map dual_cm;
	// Compute the dual of the Combinatorial map cm and return a handle to the infinite cell
	Dart_handle dh_dual = cm.dual(dual_cm, dh);
	//std::cout << "After dual function...\n";
	//display_finite_characteristics(dual_cm, DT, cell_to_dart);
	//dual_cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << dual_cm.is_without_boundary() << ", valid=" << dual_cm.is_valid() << std::endl;;
	//Update the darts and vertices of Voronoi cells
	transform_dart_to_their_dual(cm, dual_cm, cell_to_dart);
	set_geometry_of_dual(dual_cm, DT, cell_to_dart);
	//std::cout << "After setting the geometry...\n";
	//display_finite_characteristics(dual_cm, DT, cell_to_dart);
	//dual_cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << dual_cm.is_without_boundary() << ", valid=" << dual_cm.is_valid() << std::endl;;
	//std::cout << "\nDual geometry set...\n";
	removeInfinteCellsFromDual(dual_cm, dh_dual);
	//std::cout << "After Infinite cell removal...\n";
	//display_finite_characteristics(dual_cm, DT, cell_to_dart);
	//std::cout << "Finite Voronoi Diagram of " << dl_points.size() << " vertices without infinite vertex" << std::endl;
	//dual_cm.display_characteristics(std::cout) << " #withoutBoundary=AllNonFree=" << dual_cm.is_without_boundary() << ", valid=" << dual_cm.is_valid() << std::endl;;

	return dual_cm;
}

void GeometryStructure::transformCombinatorialMap(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, CGAL::Aff_transformation_3<Kernel> &transform_matrix) {

	// Extract the vertices by marking the starting edge for every 2-cell
	// Iterate every vertex (0-cell)

	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<0>::iterator face_dart = lc_combinatorial_map.one_dart_per_cell<0>().begin(),
		face_dart_end = lc_combinatorial_map.one_dart_per_cell<0>().end(); face_dart != face_dart_end; ++face_dart)
	{
		Point_CM p = lc_combinatorial_map.point(face_dart);
		//std::cout << "Point: " << p << std::endl;
		Kernel::Point_3 p3(p[0], p[1], p[2]);
		//std::cout << "Transform..." << std::endl;
		//std::cout << p3 << std::endl;
		Kernel::Point_3 p_transf(p3.transform(transform_matrix));
		Point_CM p_cm(CGAL::to_double(p_transf[0]), CGAL::to_double(p_transf[1]), CGAL::to_double(p_transf[2]));
		lc_combinatorial_map.set_vertex_attribute(face_dart, lc_combinatorial_map.create_vertex_attribute(p_cm));
	}
}

Nef_polyhedron_Ext GeometryStructure::getMovedPolyhedronPlanesFurther(Polyhedron &P, double offset) {

	std::transform(P.facets_begin(), P.facets_end(), P.planes_begin(), Plane_equation(offset));
	Nef_polyhedron_Ext new_P(Nef_polyhedron_Ext::COMPLETE);
	for (Polyhedron::Plane_iterator plane_it = P.planes_begin(); plane_it != P.planes_end(); ++plane_it) {
		Nef_polyhedron_Ext N(Nef_polyhedron_Ext::Plane_3(plane_it->a(), plane_it->b(), plane_it->c(), plane_it->d()), Nef_polyhedron_Ext::INCLUDED);
		new_P *= N;
	}
	return new_P;
}

Nef_polyhedron_Ext GeometryStructure::getBiggerPolyhedronByDistance(Linear_Complex_Combinatorial_Map lccm, double offset) {
	std::vector<Point_CM> points;
	std::vector<std::vector<int>> faces;
	std::vector<Kernel::Plane_3> planes;

	getPointsAndFaces(lccm, points, faces);

	for (auto& f : faces) {
		CGAL::Vector_3<Kernel> v_sum(0.f, 0.f, 0.f);
		for (auto& p_i : f) {
			v_sum += CGAL::Vector_3<Kernel>(points[p_i][0], points[p_i][1], points[p_i][2]);
		}

		v_sum /= CGAL::to_double(f.size());
		CGAL::Vector_3<Kernel> v_sum_unit(v_sum / sqrt(CGAL::to_double(v_sum.squared_length()))),
			v_offset(v_sum_unit * offset),
			v_parallel_plane(v_sum + v_offset);
		Kernel::Point_3 p_parallel_plane(v_parallel_plane[0], v_parallel_plane[1], v_parallel_plane[2]);
		Kernel::Plane_3 new_plane(p_parallel_plane, v_sum);
		planes.push_back(new_plane);
	}

	Nef_polyhedron_Ext new_P(Nef_polyhedron_Ext::COMPLETE);
	for (auto& plane : planes) {
		Nef_polyhedron_Ext N(Nef_polyhedron_Ext::Plane_3(plane.a(), plane.b(), plane.c(), plane.d()), Nef_polyhedron_Ext::INCLUDED);
		new_P *= N;
	}
	return new_P;
}

std::vector<std::vector<int>> removeDuplicatedPointsInIndexedFaceList(std::vector<Point_CM> points, std::vector<std::vector<int>> faces) {

	std::vector<std::vector<int>> new_faces;

	for (int i = 0; i < faces.size(); ++i) {
		std::vector<int> face = faces[i], new_face_points;
		int n_face_points = face.size();
		std::vector<bool> toskip_face_points(n_face_points, false);

		for (int i = 0; i < n_face_points - 1; ++i) {
			for (int j = i + 1; j < n_face_points; ++j) {
				if (sqrt(CGAL::to_double(CGAL::squared_distance(points[face[i]], points[face[j]]))) < DUAL_DIST_ERROR) {
					toskip_face_points[j] = true;
				}
			}
		}
		for (int k = 0; k < n_face_points; ++k) {
			if (!toskip_face_points[k]) new_face_points.push_back(face[k]);
		}
		if (new_face_points.size() > 2) {
			new_faces.push_back(new_face_points);
		}
	}
	return new_faces;
}

void GeometryStructure::getPointsAndFaces(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<Point_CM>& empty_points, std::vector<std::vector<int>>& empty_faces) {
	Linear_Complex_Combinatorial_Map::size_type mark_vertex = lc_combinatorial_map.get_new_mark();
	// Point_3 to points index vector
	std::map<Point_CM, int> points_map;
	std::map<Point_CM, int>::iterator map_it;
	int curr_point_index = 0;
	std::vector<Point_CM> points;
	std::vector<std::vector<int>> faces, final_faces;

	// Extract the vertices by marking the starting edge for every 2-cell
	// Iterate every face (2-cell): one dart for every face
	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<2>::iterator face_dart = lc_combinatorial_map.one_dart_per_cell<2>().begin(),
		face_dart_end = lc_combinatorial_map.one_dart_per_cell<2>().end(); face_dart != face_dart_end; ++face_dart)
	{
		std::vector<int> curr_face;
		Dart_handle dh_edges = face_dart;
		// Iterate every edge for every 2-cell
		while (!lc_combinatorial_map.is_marked(dh_edges, mark_vertex)) {
			Point_CM p = lc_combinatorial_map.point(dh_edges);
			map_it = points_map.find(p);
			if (map_it != points_map.end()) {
				curr_face.push_back(points_map[p]);
			}
			else {
				points.push_back(Point_CM(p.x(), p.y(), p.z()));
				points_map[p] = curr_point_index;
				curr_face.push_back(curr_point_index);
				++curr_point_index;
			}

			if (dh_edges == face_dart) {
				CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_combinatorial_map, dh_edges, mark_vertex);
			}
			dh_edges = lc_combinatorial_map.beta(dh_edges, 1);
		}
		CGAL::unmark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_combinatorial_map, dh_edges, mark_vertex);
		faces.push_back(curr_face);
	}

	int n_points = points.size();
	std::vector<bool> toskip_points(n_points, false);
	std::map<int, int> doublepoint_to_point_map;
	std::map<int, int>::iterator doublepoint_to_point_map_it;
	std::vector<Point_CM> filtered_points;

	for (int i = 0; i < n_points; ++i) {
		doublepoint_to_point_map[i] = i;
	}

	for (int i = 0; i < n_points - 1; ++i) {
		if (!toskip_points[i]) {
			for (int j = i + 1; j < n_points; ++j) {
				if (!toskip_points[j]) {
					if (sqrt(CGAL::to_double(CGAL::squared_distance(points[i], points[j]))) < DUAL_DIST_ERROR) {
						toskip_points[j] = true;
						doublepoint_to_point_map[j] = i;
					}
				}
			}
		}
	}
	int new_filtered_index = 0;
	for (int k = 0; k < n_points; ++k) {
		if (!toskip_points[k]) {
			empty_points.push_back(points[k]);
			doublepoint_to_point_map[k] = new_filtered_index;
			++new_filtered_index;
		}
		else {
			doublepoint_to_point_map[k] = doublepoint_to_point_map[doublepoint_to_point_map[k]];
		}
	}

	//Iterate all faces and remove incorrect faces and points
	/*
	std::cout << "No Filtered points...\n";
	for (auto& p : points) {
		std::cout << p << std::endl;
	}
	*/
	final_faces = removeDuplicatedPointsInIndexedFaceList(points, faces);
	/*
	std::cout << "Filtered Faces, not modified\n";

	int j = 1;
	for (auto& f : final_faces) {
		std::cout << "Face " << j << std::endl;
		for (int i = 0; i < f.size(); ++i) {
			std::cout << f[i] << " ";
		}
		std::cout << std::endl;
		++j;
	}
	*/
	for (int i = 0; i < final_faces.size(); ++i) {
		for (int j = 0; j < final_faces[i].size(); ++j) {
			final_faces[i][j] = doublepoint_to_point_map[final_faces[i][j]];
		}
	}
	/*
	std::cout << "Filtered and modified structs...\n";
	for (auto& p : empty_points) {
		std::cout << p << std::endl;
	}

	j = 1;
	for (auto& f : final_faces) {
		std::cout << "Face " << j << std::endl;
		for (int i = 0; i < f.size(); ++i) {
			std::cout << f[i] << " ";
		}
		std::cout << std::endl;
		++j;
	}
	*/
	empty_faces.assign(final_faces.begin(), final_faces.end());
}

void GeometryStructure::getPolyhedraData(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<std::vector<Point_CM>>& polyhedra_points, std::vector<std::vector<std::vector<int>>>& polyhedra_faces) {
	Linear_Complex_Combinatorial_Map::size_type mark_vertex = lc_combinatorial_map.get_new_mark();
	Linear_Complex_Combinatorial_Map::size_type mark_face = lc_combinatorial_map.get_new_mark();
	
	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<3>::iterator volume_dart = lc_combinatorial_map.one_dart_per_cell<3>().begin(),
		volume_dart_end = lc_combinatorial_map.one_dart_per_cell<3>().end(); volume_dart != volume_dart_end; ++volume_dart)
	{
		std::vector<Point_CM> points, final_points;
		std::vector<std::vector<int>> faces, final_faces;
		// Point_3 to points index vector
		std::map<Point_CM, int> points_map;
		std::map<Point_CM, int>::iterator map_it;
		int curr_point_index = 0;
		// Extract the vertices by marking the starting edge for every 2-cell
		// Iterate all darts (1-cell) checking if the face has been visited
		for (Linear_Complex_Combinatorial_Map::Dart_of_cell_range<3>::iterator face_dart = lc_combinatorial_map.darts_of_cell<3>(volume_dart).begin(),
			face_dart_end = lc_combinatorial_map.darts_of_cell<3>(volume_dart).end(); face_dart != face_dart_end; ++face_dart)
		{
			if (lc_combinatorial_map.is_marked(face_dart, mark_face)) continue;

			std::vector<int> curr_face;
			Dart_handle dh_edges = face_dart;
			bool inf_face = false;

			// Iterate every edge for every 2-cell
			while (!lc_combinatorial_map.is_marked(dh_edges, mark_vertex)) {
				Point_CM p = lc_combinatorial_map.point(dh_edges);
				map_it = points_map.find(p);
				if (map_it != points_map.end()) {
					curr_face.push_back(points_map[p]);
				}
				else {
					if ( (CGAL::is_finite(p[0]) && (CGAL::is_finite(p[1])) && (CGAL::is_finite(p[2]))) ) {
						points.push_back(Point_CM(p[0], p[1], p[2]));
						points_map[p] = curr_point_index;
						curr_face.push_back(curr_point_index);
						++curr_point_index;
					}
					// A point of a face is infinite
					else {
						inf_face = true;
					}
				}
				if (dh_edges == face_dart) {
					lc_combinatorial_map.mark(dh_edges, mark_vertex);
				}
				lc_combinatorial_map.mark(dh_edges, mark_face);

				dh_edges = lc_combinatorial_map.beta(dh_edges, 1);
			}
			lc_combinatorial_map.unmark(dh_edges, mark_vertex);
			if (!inf_face)
				faces.push_back(curr_face);
		}		

		for (Linear_Complex_Combinatorial_Map::Dart_of_cell_range<3>::iterator face_dart = lc_combinatorial_map.darts_of_cell<3>(volume_dart).begin(),
			face_dart_end = lc_combinatorial_map.darts_of_cell<3>(volume_dart).end(); face_dart != face_dart_end; ++face_dart)
		{
			lc_combinatorial_map.unmark(face_dart, mark_face);
		}

		int n_points = points.size();
		std::vector<bool> toskip_points(n_points, false);
		std::map<int, int> doublepoint_to_point_map;
		std::vector<Point_CM> filtered_points;

		for (int i = 0; i < n_points; ++i) {
			doublepoint_to_point_map[i] = i;
		}

		for (int i = 0; i < n_points - 1; ++i) {
			if (!toskip_points[i]) {
				for (int j = i + 1; j < n_points; ++j) {
					if (!toskip_points[j]) {
						if (sqrt(CGAL::to_double(CGAL::squared_distance(points[i], points[j]))) < DUAL_DIST_ERROR) {
							toskip_points[j] = true;
							doublepoint_to_point_map[j] = i;
						}
					}
				}
			}
		}
		int new_filtered_index = 0;
		for (int k = 0; k < n_points; ++k) {
			if (!toskip_points[k]) {
				final_points.push_back(points[k]);
				doublepoint_to_point_map[k] = new_filtered_index;
				++new_filtered_index;
			}
			else {
				doublepoint_to_point_map[k] = doublepoint_to_point_map[doublepoint_to_point_map[k]];
			}
		}

		//Iterate all faces and remove incorrect faces and points
		final_faces = removeDuplicatedPointsInIndexedFaceList(points, faces);
				
		for (int i = 0; i < final_faces.size(); ++i) {
			for (int j = 0; j < final_faces[i].size(); ++j) {
				final_faces[i][j] = doublepoint_to_point_map[final_faces[i][j]];
			}
		}
				
		polyhedra_faces.push_back(final_faces);
		polyhedra_points.push_back(final_points);
	}

}

void GeometryStructure::getPolyhedraData_New(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<std::vector<Point_CM>>& polyhedra_points, std::vector<std::vector<std::vector<int>>>& polyhedra_faces, Eigen::VectorXd& limit_inf, Eigen::VectorXd& limit_sup) {
	Linear_Complex_Combinatorial_Map::size_type mark_vertex = lc_combinatorial_map.get_new_mark();
	Linear_Complex_Combinatorial_Map::size_type mark_face = lc_combinatorial_map.get_new_mark();

	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<3>::iterator volume_dart = lc_combinatorial_map.one_dart_per_cell<3>().begin(),
		volume_dart_end = lc_combinatorial_map.one_dart_per_cell<3>().end(); volume_dart != volume_dart_end; ++volume_dart)
	{
		std::vector<Point_CM> points, final_points;
		std::vector<std::vector<int>> faces, final_faces;
		// Point_3 to points index vector
		std::map<Point_CM, int> points_map;
		std::map<Point_CM, int>::iterator map_it;
		int curr_point_index = 0, point_counter = 0;
		bool inf_poly = false;

		// Extract the vertices by marking the starting edge for every 2-cell
		// Iterate all darts (1-cell) checking if the face has been visited
		for (Linear_Complex_Combinatorial_Map::Dart_of_cell_range<3>::iterator face_dart = lc_combinatorial_map.darts_of_cell<3>(volume_dart).begin(),
			face_dart_end = lc_combinatorial_map.darts_of_cell<3>(volume_dart).end(); face_dart != face_dart_end; ++face_dart)
		{
			if (lc_combinatorial_map.is_marked(face_dart, mark_face)) continue;

			std::vector<int> curr_face;
			Dart_handle dh_edges = face_dart;

			// Iterate every edge for every 2-cell
			while (!lc_combinatorial_map.is_marked(dh_edges, mark_vertex)) {
				Point_CM p = lc_combinatorial_map.point(dh_edges);
				map_it = points_map.find(p);
				if (map_it != points_map.end()) {
					curr_face.push_back(points_map[p]);
				}
				else if ((CGAL::is_finite(p[0]) && p[0] >= limit_inf[0] - DOMAIN_LIMIT && p[0] <= limit_sup[0] + DOMAIN_LIMIT) &&
						(CGAL::is_finite(p[1]) && p[1] >= limit_inf[1] - DOMAIN_LIMIT && p[1] <= limit_sup[1] + DOMAIN_LIMIT) &&
						(CGAL::is_finite(p[2]) && p[2] >= limit_inf[2] - DOMAIN_LIMIT && p[2] <= limit_sup[2] + DOMAIN_LIMIT)) {
					points.push_back(Point_CM(p[0], p[1], p[2]));
					points_map[p] = curr_point_index;
					curr_face.push_back(curr_point_index);
					++curr_point_index;
					++point_counter;
				}
				// A point of a face is infinite or out of boundary
				else {
					inf_poly = true;
				}
				if (dh_edges == face_dart) {
					lc_combinatorial_map.mark(dh_edges, mark_vertex);
				}
				lc_combinatorial_map.mark(dh_edges, mark_face);

				dh_edges = lc_combinatorial_map.beta(dh_edges, 1);
			}
			lc_combinatorial_map.unmark(dh_edges, mark_vertex);
			if (!inf_poly)
				faces.push_back(curr_face);
			else {
				curr_point_index -= point_counter;
				break;
			}
		}

		for (Linear_Complex_Combinatorial_Map::Dart_of_cell_range<3>::iterator face_dart = lc_combinatorial_map.darts_of_cell<3>(volume_dart).begin(),
			face_dart_end = lc_combinatorial_map.darts_of_cell<3>(volume_dart).end(); face_dart != face_dart_end; ++face_dart)
		{
			lc_combinatorial_map.unmark(face_dart, mark_face);
		}

		if (inf_poly) {
			continue;
		}

		int n_points = points.size();
		std::vector<bool> toskip_points(n_points, false);
		std::map<int, int> doublepoint_to_point_map;
		std::vector<Point_CM> filtered_points;

		for (int i = 0; i < n_points; ++i) {
			doublepoint_to_point_map[i] = i;
		}

		for (int i = 0; i < n_points - 1; ++i) {
			if (!toskip_points[i]) {
				for (int j = i + 1; j < n_points; ++j) {
					if (!toskip_points[j]) {
						if (sqrt(CGAL::to_double(CGAL::squared_distance(points[i], points[j]))) < DUAL_DIST_ERROR) {
							toskip_points[j] = true;
							doublepoint_to_point_map[j] = i;
						}
					}
				}
			}
		}
		int new_filtered_index = 0;
		for (int k = 0; k < n_points; ++k) {
			if (!toskip_points[k]) {
				final_points.push_back(points[k]);
				doublepoint_to_point_map[k] = new_filtered_index;
				++new_filtered_index;
			}
			else {
				doublepoint_to_point_map[k] = doublepoint_to_point_map[doublepoint_to_point_map[k]];
			}
		}

		//Iterate all faces and remove incorrect faces and points
		final_faces = removeDuplicatedPointsInIndexedFaceList(points, faces);

		for (int i = 0; i < final_faces.size(); ++i) {
			for (int j = 0; j < final_faces[i].size(); ++j) {
				final_faces[i][j] = doublepoint_to_point_map[final_faces[i][j]];
			}
		}

		polyhedra_faces.push_back(final_faces);
		polyhedra_points.push_back(final_points);
	}

}

Polyhedron GeometryStructure::getPolyhedronFromCombinatorialMap(Linear_Complex_Combinatorial_Map lc_combinatorial_map, bool subdivide) {
	Polyhedron P;
	//Triangulation T;
	//Triangulation_basic T;
	//Build_Polyhedron_From_Triangulation<Polyhedron::HalfedgeDS> poly(T);
	Build_Polyhedron_From_Combinatorial_Map<Polyhedron::HalfedgeDS> poly(lc_combinatorial_map);
	P.delegate(poly);
	P.normalize_border();
	if (subdivide) {
		if (P.size_of_border_edges() != 0) {
			std::cerr << "The input object has border edges. Cannot subdivide."
				<< std::endl;
			std::exit(1);
		}
		subdiv(P);
	}

	return P;
}

void GeometryStructure::writePointsToCSVFormatFile(const char* filename) {

	std::ofstream out_points;
	out_points.open(filename);
	out_points << "X,Y,Z" << std::endl;
	for (Eigen::VectorXd &p : this->getPoints()) {
		for (int i = 0; i < 3; ++i) {
			if (i == 2) {
				out_points << p(i);
			}
			else {
				out_points << p(i) << ",";
			}
		}
		out_points << std::endl;
	}
	out_points.close();
}
