/*Copyright (C) 2018-2022,2026 Marco M. Mosca

This file is part of VoronoiLatticeDistances.

VoronoiLatticeDistances is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VoronoiLatticeDistances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VoronoiLatticeDistances. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef _LINEARCOMPLEX_H
#define _LINEARCOMPLEX_H
#define INF_COMPONENT 1e+4
#define VOR_INF_COMPONENT 1e+2
#define DOMAIN_LIMIT 20
#define DUAL_DIST_ERROR 0.05

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <geom.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> Linear_Complex_Combinatorial_Map;
typedef Linear_Complex_Combinatorial_Map::Dart_handle	Dart_handle;
typedef Linear_Complex_Combinatorial_Map::Point	Point_CM;
typedef CGAL::Delaunay_triangulation_3<Linear_Complex_Combinatorial_Map::Traits> Delaunay_Triangulation;
typedef CGAL::Triangulation_3<Linear_Complex_Combinatorial_Map::Traits> Triangulation;
typedef CGAL::Triangulation_3<Kernel> Triangulation_basic;

typedef Delaunay_Triangulation::All_cells_iterator  ACell_iterator;
typedef Delaunay_Triangulation::All_facets_iterator  AFacet_iterator;
typedef Delaunay_Triangulation::Finite_cells_iterator FCell_iterator;
typedef Delaunay_Triangulation::Finite_facets_iterator FFacet_iterator;
typedef Triangulation::All_cells_iterator T_ACell_iterator;
typedef Triangulation::All_facets_iterator T_AFacet_iterator;
typedef Triangulation::Finite_facets_iterator T_FFacet_iterator;
typedef Triangulation::All_vertices_iterator T_AVertex_iterator;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation_basic::All_cells_iterator Tb_ACell_iterator;
typedef Triangulation_basic::All_facets_iterator Tb_AFacets_iterator;
typedef Triangulation_basic::Cell_handle Tb_Cell_Handle;




void transform_dart_to_their_dual(Linear_Complex_Combinatorial_Map& lc_combinatorial_map, 
	Linear_Complex_Combinatorial_Map& lc_dual_combinatorial_map, 
	std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart);

void set_geometry_of_dual(Linear_Complex_Combinatorial_Map& lc_dual_combinatorial_map,
	Delaunay_Triangulation& lc_combinatorial_map,
	std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart);

void removeInfinteCellsFromDual(Linear_Complex_Combinatorial_Map &dual_cm, Dart_handle dh_inf);

/*
*	Remove all infinite simplices from a Linear Complex Combinatorial Map
*/
void updateToFiniteSimplicialComplex(Linear_Complex_Combinatorial_Map& lc_combinatorial_map,
	Delaunay_Triangulation& DT,
	std::map<typename Delaunay_Triangulation::Cell_handle, typename Linear_Complex_Combinatorial_Map::Dart_handle>& cell_to_dart);


/*
*	Remove all faces which have incorrect geometry from a Linear Complex Combinatorial Map
*/
void removeIncorrectGeometryFaces(Linear_Complex_Combinatorial_Map& lc_combinatorial_map);

#endif // !_LINEARCOMPLEX_H
