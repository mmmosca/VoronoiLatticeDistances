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
#ifndef _KERNEL_H
#define _KERNEL_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<double> Kernel_simple;
typedef CGAL::Extended_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Kernel_Ext;

#endif
