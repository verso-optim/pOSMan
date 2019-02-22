/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <cassert>
#include <iostream>
#include <numeric>
#include <unistd.h>

#include "utils/osm_parser.h"

void display_usage() {
  std::string usage = "Copyright (C) 2019, VERSO\n";
  usage += "\n\tposman [OPTION]... -i FILE\n";
  usage += "Options:\n";
  usage += "\t-e EDGES,\t\t\t file containing OSM edges\n";
  usage += "\t-n NODES,\t\t\t file containing OSM nodes\n";
  std::cout << usage << std::endl;
  exit(0);
}

int main(int argc, char** argv) {
  // Parsing command-line arguments.
  const char* optString = "e:n:h?";
  int opt = getopt(argc, argv, optString);

  std::string edges_file;
  std::string nodes_file;

  while (opt != -1) {
    switch (opt) {
    case 'e':
      edges_file = optarg;
      break;
    case 'h':
      display_usage();
      break;
    case 'n':
      nodes_file = optarg;
      break;
    default:
      break;
    }
    opt = getopt(argc, argv, optString);
  }

  std::cout << "[info] Loading graph." << std::endl;
  auto graph = posman::io::parse(nodes_file, edges_file);

  unsigned adjacencies =
    std::accumulate(graph.adjacency_list.begin(),
                    graph.adjacency_list.end(),
                    0,
                    [&](auto sum, const auto& l) { return sum + l.size(); });
  assert(adjacencies % 2 == 0);

  std::cout << "  Graph has " << graph.nodes.size() << " nodes and "
            << (adjacencies / 2) << " edges." << std::endl;

  return 0;
}
