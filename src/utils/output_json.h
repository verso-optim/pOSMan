#ifndef OUTPUT_JSON_H
#define OUTPUT_JSON_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <string>

#include "../include/rapidjson/document.h"

#include "structures/graph/undirected_graph.h"

namespace posman {
namespace io {

void log_graph_as_geojson(const UndirectedGraph& graph,
                          const std::string& output_file);

void write_output(const UndirectedGraph& graph,
                  const std::vector<Index>& path,
                  const std::vector<Id>& path_way_ids,
                  const std::string& output_file);

} // namespace io
} // namespace posman

#endif
