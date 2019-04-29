#ifndef OSM_PARSER_H
#define OSM_PARSER_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include "structures/graph/undirected_graph.h"

namespace posman {
namespace io {

UndirectedGraph parse_graph(const std::string& nodes_filename,
                            const std::string& edges_filename,
                            GeometryList& geometries);

UndirectedGraph parse_ways(const std::string& ways_filename,
                           const UndirectedGraph& global_graph);

} // namespace io
} // namespace posman

#endif
