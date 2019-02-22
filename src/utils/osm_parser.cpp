/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <cassert>
#include <fstream>
#include <unordered_map>

#include "structures/graph/node.h"
#include "utils/osm_parser.h"

namespace posman {
namespace io {

UndirectedGraph parse(const std::string& nodes_filename,
                      const std::string& edges_filename) {
  UndirectedGraph graph;

  // Parsing OSM nodes.
  std::unordered_map<Id, Node> nodes;

  std::ifstream nodes_file(nodes_filename);

  std::string current_line;
  std::getline(nodes_file, current_line); // Unused header line.

  while (std::getline(nodes_file, current_line)) {
    Node current_node(current_line);
    nodes.insert(std::make_pair(current_node.osm_id, current_node));
  }

  // Parsing OSM ways.
  std::ifstream edges_file(edges_filename);
  std::getline(edges_file, current_line); // Unused header line.

  while (std::getline(edges_file, current_line)) {
    std::string::size_type e = current_line.find(',');
    Id osm_way_id = strtoul(current_line.substr(0, e).c_str(), nullptr, 10);

    std::string::size_type s = e + 1;
    e = current_line.find(',', s);
    Id source_node_id =
      strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);
    s = e + 1;
    e = current_line.find(',', s);
    Id target_node_id =
      strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    s = e + 1;
    e = current_line.find(',', s);
    Distance length =
      static_cast<Distance>(100 * std::stod(current_line.substr(s, e - s)));

    // Profile statuses.
    s = e + 1;
    e = current_line.find(',', s);
    auto foot = strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);
    // // car_forward
    // s = e + 1;
    // e = current_line.find(',', s);
    // auto car_forward = strtoul(current_line.substr(s, e - s).c_str(),
    // nullptr, 10);
    // // car_backward
    // s = e + 1;
    // e = current_line.find(',', s);
    // auto car_backward = strtoul(current_line.substr(s, e - s).c_str(),
    // nullptr, 10);
    // // bike_forward
    // s = e + 1;
    // e = current_line.find(',', s);
    // auto bike_forward = strtoul(current_line.substr(s, e - s).c_str(),
    // nullptr, 10);
    // // bike_backward
    // s = e + 1;
    // e = current_line.find(',', s);
    // auto bike_backward = strtoul(current_line.substr(s, e - s).c_str(),
    // nullptr, 10);

    if (foot != 0) {
      auto source = nodes.find(source_node_id);
      assert(source != nodes.end());

      auto target = nodes.find(target_node_id);
      assert(target != nodes.end());

      graph.add_edge(osm_way_id, source->second, target->second, length);
    }
  }

  return graph;
}

} // namespace io
} // namespace posman
