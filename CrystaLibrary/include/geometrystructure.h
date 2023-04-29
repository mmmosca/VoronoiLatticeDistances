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
#ifndef _GEOMETRYSTRUCTURE_H
#define _GEOMETRYSTRUCTURE_H

#include <Eigen/Dense>
#include <vector>
#include <periodicpointcloud.h>
#include <linearcomplex.h>
#include <kernel.h>

template<class HDS>
class Build_Polyhedron_From_Triangulation : public CGAL::Modifier_base<HDS> {
public:
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point_3 Point_3;
	std::vector<Kernel::Point_3> poly_points;
	std::vector<std::vector<int>> poly_faces;
	Triangulation poly_triangulation;
	Build_Polyhedron_From_Triangulation(Triangulation &T) {
		poly_triangulation = T;
	}
	void operator()(HDS& hds) {
		CGAL::Polyhedron_incremental_builder_3<HDS> Poly_Builder(hds, true);

		// Convert the triangulation to a Linear Complex Combinatorial Map
		Linear_Complex_Combinatorial_Map cm;
		// Make a map from the cell (3-simplex : volume) -> dart handle (all the volumes have a dart handle)
		std::map<Triangulation::Cell_handle, Linear_Complex_Combinatorial_Map::Dart_handle> cell_to_dart;
		// Update the combinatorial map: Tetrahedra are added to the Linear Complex combinatorial map (No change in darts). Returns a dart handle to an infinite cell (3-simplex)
		Dart_handle dh = CGAL::import_from_triangulation_3<Linear_Complex_Combinatorial_Map, Triangulation>(cm, poly_triangulation, &cell_to_dart);
		//Linear_Complex_Combinatorial_Map::size_type mark_vertex = cm.get_new_mark();		
		removeIncorrectGeometryFaces(cm);
		// Collect points and indexes before calling Poly_builder.begin_surface()
		// Add check on infinite cell
		GeometryStructure::getPointsAndFaces(cm, poly_points, poly_faces);
		Poly_Builder.begin_surface(poly_points.size(), poly_faces.size());
		for (int j = 0; j < poly_points.size(); ++j) {
			Point_3 p = Point_3(poly_points[j][0], poly_points[j][1], poly_points[j][2]);
			Poly_Builder.add_vertex(p);
		}
		std::vector<int> curr_face;
		int face_size;
		for (int i = 0; i < poly_faces.size(); ++i) {
			curr_face = poly_faces[i];
			face_size = curr_face.size();

			Poly_Builder.begin_facet();
			for (int j = 0; j < face_size; ++j) {
				Poly_Builder.add_vertex_to_facet(curr_face[j]);
			}
			Poly_Builder.end_facet();
		}
		Poly_Builder.end_surface();

	}
};

template<class HDS>
class Build_Polyhedron_From_Combinatorial_Map : public CGAL::Modifier_base<HDS> {
public:
	typedef typename HDS::Vertex   Vertex;
	typedef typename Vertex::Point_3 Point_3;
	std::vector<Point_CM> poly_points;
	std::vector<std::vector<int>> poly_faces;
	Linear_Complex_Combinatorial_Map lcc;
	Build_Polyhedron_From_Combinatorial_Map(Linear_Complex_Combinatorial_Map& lc_cm) {
		lcc = lc_cm;
	}
	void operator()(HDS& hds) {
		CGAL::Polyhedron_incremental_builder_3<HDS> Poly_Builder(hds, true);
		// Collect points and indexes before calling Poly_builder.begin_surface()
		// Add check on infinite cell
		GeometryStructure::getPointsAndFaces(lcc, poly_points, poly_faces);
		int halfedges_number = 0;
		for (auto& f : poly_faces) {
			halfedges_number += f.size();
		}
		Poly_Builder.begin_surface(poly_points.size(), poly_faces.size(), 0, Poly_Builder.ABSOLUTE_INDEXING);
		for (int j = 0; j < poly_points.size(); ++j) {
			Point_3 p = Point_3(poly_points[j][0], poly_points[j][1], poly_points[j][2]);
			Poly_Builder.add_vertex(p);
		}
		std::vector<int> curr_face;
		int face_size;
		for (int i = 0; i < poly_faces.size(); ++i) {
			curr_face = poly_faces[i];
			face_size = curr_face.size();

			Poly_Builder.begin_facet();
			for (int j = 0; j < face_size; ++j) {
				Poly_Builder.add_vertex_to_facet(curr_face[j]);
			}
			Poly_Builder.end_facet();
		}
		Poly_Builder.remove_unconnected_vertices();
		Poly_Builder.end_surface();
	}
};

struct Plane_equation {
	double offset;
	Plane_equation(double off = 0) : offset(off) {};
	Polyhedron::Facet::Plane_3 operator()(Polyhedron::Facet& f) {
		typename Polyhedron::Facet::Halfedge_handle h = f.halfedge();

		// Circulators for vertices
		Polyhedron::Halfedge_around_facet_circulator v_circ;
		v_circ = f.facet_begin();
		CGAL::Vector_3<Kernel> v_sum(0.f, 0.f, 0.f);
		do {
			Kernel::Point_3 p = v_circ->vertex()->point();
			v_sum += CGAL::Vector_3<Kernel>(p[0], p[1], p[2]);
		} while (++v_circ != f.facet_begin());

		v_sum /= CGAL::to_double(f.size());

		CGAL::Vector_3<Kernel> v_sum_unit(v_sum / sqrt(CGAL::to_double(v_sum.squared_length()))),
			v_offset(v_sum_unit * offset),
			v_parallel_plane(v_sum + v_offset);
		Kernel::Point_3 p_parallel_plane(v_parallel_plane[0], v_parallel_plane[1], v_parallel_plane[2]);
		Kernel::Plane_3 new_plane(p_parallel_plane, v_sum);
		//std::cout << "V_sum: " << v_sum << " V_sum_unit: " << v_sum_unit << " V_offset: " << v_offset  << " offset:" << sqrt(CGAL::to_double(v_offset.squared_length())) << "=" << offset << std::endl;


		return new_plane;
	}
};

class GeometryStructure : public PeriodicPointCloud {
public:
	GeometryStructure() : PeriodicPointCloud(3) {};

	std::vector<Eigen::VectorXd> getPoints();

	// Get The Delaunay Triangulation with infinite vertex of the spanned lattice
	Linear_Complex_Combinatorial_Map getDelaunayTriangulation();

	// Get The finite Delaunay Triangulation of the spanned lattice
	Linear_Complex_Combinatorial_Map getFiniteDelaunayTriangulation();

	// Get The Voronoi diagram with infinite verteces and the handle to an infinite cell
	Linear_Complex_Combinatorial_Map getVoronoiDiagram(Dart_handle &dh_inf);

	// Get The finite Voronoi diagram of the spanned lattice
	Linear_Complex_Combinatorial_Map getFiniteVoronoiDiagram();

	static void transformCombinatorialMap(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, CGAL::Aff_transformation_3<Kernel>& transform_matrix);

	static Nef_polyhedron_Ext getMovedPolyhedronPlanesFurther(Polyhedron &P, double offset);

	static Nef_polyhedron_Ext getBiggerPolyhedronByDistance(Linear_Complex_Combinatorial_Map lccm, double offset);

	static void getPointsAndFaces(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<Point_CM>& empty_points, std::vector<std::vector<int>>& empty_faces);

	static void getPolyhedraData(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<std::vector<Point_CM>>& polyhedra_points, std::vector<std::vector<std::vector<int>>>& polyhedra_faces);

	static void getPolyhedraData_New(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, std::vector<std::vector<Point_CM>>& polyhedra_points, std::vector<std::vector<std::vector<int>>>& polyhedra_faces, Eigen::VectorXd& limit_inf, Eigen::VectorXd& limit_sup);

	static Polyhedron getPolyhedronFromCombinatorialMap(Linear_Complex_Combinatorial_Map lc_combinatorial_map, bool subdivide = true);

	/*
	*	It is required that the Combinatorial maps have their centre at the origin (0,0,0)
	*/
	void writePointsToCSVFormatFile(const char* filename);
};


#endif //!_GEOMETRYSTRUCTURE_H