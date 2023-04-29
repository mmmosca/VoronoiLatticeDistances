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
#ifndef _GEOM_H
#define _GEOM_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <assert.h>
#include <queue>
#include <map>
#include <kernel.h>

#include <boost/math/special_functions.hpp>
#include <boost/math/constants/constants.hpp>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Side_of_triangle_mesh.h>


/*Used to define basic geometric objects*/
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel_Ext>  Nef_polyhedron_Ext;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef CGAL::Polyhedron_3<Kernel_Ext>  Polyhedron_Ext;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
typedef CGAL::Polytope_distance_d_traits_3<Kernel> Polytope_distance_d_Traits;
typedef CGAL::Polytope_distance_d<Polytope_distance_d_Traits>     Polytope_distance;
typedef CGAL::Side_of_triangle_mesh<CGAL::Polyhedron_3<Kernel>, Kernel> Point_inside_Poly3;


#define M_PI           3.14159265358979323846
#define ROUND_D	20
#define PRINTSTEPINFO(step, r1, v1, r2, v2, r3, v3)	std::cout << "|--> Step " << (step) << ":\n" << \
													"\t" << to_string((r1)) << ":\t" << (v1)->transpose() << ",\n" << \
													"\t" << to_string((r2)) << ":\t" << (v2)->transpose() << ",\n" << \
													"\t" << to_string((r3)) << ":\t" << (v3)->transpose() << std::endl;
#define PRINTREDUCEINFO(message, r1, r2, k)	std::cout << "\t|--> " << message << ": " << (r1) << "->" << (r2) << ", " << (k) << std::endl;
#define math_sign(X) ((X) < 0.f) ? -1 : 1
#define MAX(X,Y)	((X)>=(Y)) ? (X) : (Y)
#define MIN(X,Y)	((X)<=(Y)) ? (X) : (Y)

enum class Axis {X,Y,Z};

using namespace std;

double roundToNthDecimal(double i, int n);

/* 
*	This function returns the Cartesian vector components of the 3 unit cell axis 
*	param1 : vector of lengths (a, b, c) and angles (alpha, beta, gamma)
*/
Eigen::Matrix3d getTransformationMatrixFromFractionalToCartesian(std::vector<double> params);

double getEuclideanDistance(Eigen::VectorXd p1, Eigen::VectorXd p2);

Eigen::Matrix3d getRotationMatrix(Axis a, double angle_degree);

Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d a, double angle_rad );

Kernel::Vector_3 getRotatedVector(Kernel::Vector_3 p, Kernel::Vector_3 a, double angle_rad);

Kernel::Point_3 getRotatedPoint(Kernel::Vector_3 p, Kernel::Vector_3 a, double angle_rad);

Kernel_simple::Point_3 getRotatedPoint(Kernel_simple::Vector_3 p, Kernel_simple::Vector_3 a, double angle_rad);

Eigen::Matrix3d getRotationMatricesComposition(double theta, double phi, double psi);

void processNewFractionalCoordinates(Eigen::VectorXd fractional, Eigen::VectorXd& new_fractional);

/********************\
|**** POLYHEDRON ****|
\********************/

struct Smooth_old_vertex {
	Kernel::Point_3 operator()(const Polyhedron::Vertex& v) const {
		CGAL_precondition((CGAL::circulator_size(v.vertex_begin()) & 1) == 0);
		std::size_t degree = CGAL::circulator_size(v.vertex_begin()) / 2;
		double alpha = (4.0 - 2.0 * std::cos(2.0 * CGAL_PI / degree)) / 9.0;
		Kernel::Vector_3 vec = (v.point() - CGAL::ORIGIN) * (1.0 - alpha);
		Polyhedron::Halfedge_around_vertex_const_circulator h = v.vertex_begin();
		do {
			vec = vec + (h->opposite()->vertex()->point() - CGAL::ORIGIN)
				* alpha / static_cast<double>(degree);
			++h;
			CGAL_assertion(h != v.vertex_begin()); // even degree guaranteed
			++h;
		} while (h != v.vertex_begin());
		return (CGAL::ORIGIN + vec);
	}
};

void create_center_vertex(Polyhedron& P, Polyhedron::Facet_iterator f);

void flip_edge(Polyhedron& P, Polyhedron::Halfedge_handle e);

void subdiv(Polyhedron& P);
#endif // !_GEOM_H
