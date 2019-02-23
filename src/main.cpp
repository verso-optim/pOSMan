/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unistd.h>

#include "structures/typedefs.h"
#include "utils/osm_parser.h"
#include "utils/output_json.h"

void display_usage() {
  std::string usage = "Copyright (C) 2019, VERSO\n";
  usage += "\n\tposman [OPTION]... -i FILE\n";
  usage += "Options:\n";
  usage += "\t-e EDGES,\t\t\t file containing OSM edges\n";
  usage += "\t-g GEOJSON,\t\t\t file to write target graph\n";
  usage += "\t-n NODES,\t\t\t file containing OSM nodes\n";
  usage += "\t-w WAYS,\t\t\t file containing target ways\n";
  std::cout << usage << std::endl;

  exit(0);
}

int main(int argc, char** argv) {
  // Parsing command-line arguments.
  const char* optString = "e:g:h?n:w:";
  int opt = getopt(argc, argv, optString);

  std::string edges_file;
  std::string geojson_target;
  std::string nodes_file;
  std::string ways_file;

  while (opt != -1) {
    switch (opt) {
    case 'e':
      edges_file = optarg;
      break;
    case 'g':
      geojson_target = optarg;
      break;
    case 'h':
      display_usage();
      break;
    case 'n':
      nodes_file = optarg;
      break;
    case 'w':
      ways_file = optarg;
      break;
    default:
      break;
    }
    opt = getopt(argc, argv, optString);
  }

  std::cout << "[info] Loading global graph." << std::endl;
  auto global_graph = posman::io::parse_graph(nodes_file, edges_file);

  std::cout << "  Global graph has " << global_graph.number_of_nodes()
            << " nodes and " << global_graph.number_of_edges() << " edges."
            << std::endl;

  std::cout << "[info] Loading target graph." << std::endl;
  auto target_graph = posman::io::parse_ways(ways_file, global_graph);

  unsigned target_nodes = target_graph.number_of_nodes();

  std::cout << "  Target graph has " << target_nodes << " nodes and "
            << target_graph.number_of_edges() << " edges." << std::endl;

  unsigned even_nodes =
    std::count_if(target_graph.adjacency_list.begin(),
                  target_graph.adjacency_list.end(),
                  [](const auto& l) { return (l.size() % 2) == 0; });

  std::cout << "  Target graph has " << even_nodes
            << " nodes with even degree ("
            << static_cast<int>(100 * static_cast<float>(even_nodes) /
                                target_nodes)
            << "%)" << std::endl;

  if (!geojson_target.empty()) {
    posman::io::log_graph_as_geojson(target_graph, geojson_target);
  }

  return 0;
}
