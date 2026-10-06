#pragma once

// Cabecera paraguas de Geometra: incluye toda la API pública.
//
//     #include <geometra/geometra.hpp>
//     auto hull = geo::convex_hull(points);
//     auto mesh = geo::delaunay_triangulation(points);

#include "geometra/convex_hull.hpp"
#include "geometra/delaunay.hpp"
#include "geometra/point.hpp"
#include "geometra/polygon.hpp"
#include "geometra/predicates.hpp"
#include "geometra/triangulation.hpp"
#include "geometra/version.hpp"
