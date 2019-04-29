/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "structures/graph/node.h"
#include "utils/osm_parser.h"

namespace posman {
namespace io {

UndirectedGraph parse_graph(const std::string& nodes_filename,
                            const std::string& edges_filename,
                            GeometryList& geometries) {
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
    // car_forward
    s = e + 1;
    e = current_line.find(',', s);
    // auto car_forward =
    //   strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    // car_backward
    s = e + 1;
    e = current_line.find(',', s);
    // auto car_backward =
    //   strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    // bike_forward
    s = e + 1;
    e = current_line.find(',', s);
    // auto bike_forward =
    //   strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    // bike_backward
    s = e + 1;
    e = current_line.find(',', s);
    // auto bike_backward =
    //   strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    // foot
    s = e + 1;
    e = current_line.find(',', s);
    auto foot = strtoul(current_line.substr(s, e - s).c_str(), nullptr, 10);

    if (foot != 0) {
      auto source = nodes.find(source_node_id);
      assert(source != nodes.end());

      auto target = nodes.find(target_node_id);
      assert(target != nodes.end());

      // Parse edge LINESTRING.
      std::vector<std::array<Coordinate, 2>> geometry;
      s = e + 1;
      assert(current_line.substr(s, 12) == "\"LINESTRING(");
      std::string::size_type end = current_line.size() - 2;
      assert(current_line.substr(end) == ")\"");

      std::istringstream linestring_stream(
        current_line.substr(s + 12, end - s - 12));
      std::string coordinates;
      while (std::getline(linestring_stream, coordinates, ',')) {
        Coordinate lon, lat;
        std::stringstream ss(coordinates);
        ss >> lon >> lat;
        geometry.push_back({lon, lat});
      }

      graph.add_edge(osm_way_id, source->second, target->second, length);
      geometries[source->second.osm_id][target->second.osm_id] =
        std::move(geometry);
    }
  }

  return graph;
}

struct custom_hash {
  template <class T1, class T2>
  std::size_t operator()(std::pair<T1, T2> const& pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

UndirectedGraph parse_ways(const std::string& ways_filename,
                           const UndirectedGraph& global_graph) {
  std::ifstream ways_file(ways_filename);
  std::string current_line;

  std::unordered_set<Id> ways_ids;
  std::unordered_set<Id> done_ways;

  while (std::getline(ways_file, current_line)) {
    ways_ids.insert(strtoul(current_line.c_str(), nullptr, 10));
  }

  UndirectedGraph target_graph;
  Distance total_length = 0;

  // Remember pair with source and target node id to avoid duplicates,
  // as edges are stored in both ways in adjacency list.
  std::unordered_set<std::pair<Id, Id>, custom_hash> source_target;

  for (std::size_t i = 0; i < global_graph.number_of_nodes(); ++i) {
    auto& source_node = global_graph.nodes[i];
    auto source_id = source_node.osm_id;

    for (const auto& edge : global_graph.adjacency_list[i]) {
      auto search_way = ways_ids.find(edge.osm_way_id);
      if (search_way == ways_ids.end()) {
        // Current edge from global graph should not be included in
        // target graph.
        continue;
      }
      done_ways.insert(edge.osm_way_id);

      auto& target_node = global_graph.nodes[edge.to];
      auto target_id = target_node.osm_id;

      auto search = source_target.find(std::make_pair(target_id, source_id));
      if (search != source_target.end()) {
        // This edge has already been added the other way around.
        continue;
      }

      target_graph.add_edge(edge.osm_way_id,
                            source_node,
                            target_node,
                            edge.length);
      total_length += edge.length;
      source_target.insert(std::make_pair(source_id, target_id));
    }
  }

  std::cout << "[Info] Total length for ways to visit: " << total_length
            << " cm." << std::endl;

  for (auto id : ways_ids) {
    if (done_ways.find(id) == done_ways.end()) {
      std::cout << "[Warning] Way id not found in OSM extract: " << id
                << std::endl;
    }
  }

  return target_graph;
}

} // namespace io
} // namespace posman
