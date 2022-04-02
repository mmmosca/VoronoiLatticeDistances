#ifndef _KERNEL_H
#define _KERNEL_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Simple_cartesian<double> Kernel_simple;
typedef CGAL::Extended_cartesian< CGAL::Lazy_exact_nt<CGAL::Gmpq> > Kernel_Ext;

#endif