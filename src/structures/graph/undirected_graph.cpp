/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <cassert>

#include "structures/graph/undirected_graph.h"

namespace posman {

bool UndirectedGraph::has_node(const Node& node) {
  return _osm_id_to_index.find(node.osm_id) != _osm_id_to_index.end();
}

void UndirectedGraph::add_node(const Node& node) {
  _osm_id_to_index.insert(std::make_pair(node.osm_id, nodes.size()));
  nodes.emplace_back(node.osm_id, node.lon, node.lat);
  adjacency_list.push_back({});
}

void UndirectedGraph::add_edge(Id osm_way_id,
                               const Node& source,
                               const Node& target,
                               Distance length) {
  if (!this->has_node(source)) {
    this->add_node(source);
  }
  if (!this->has_node(target)) {
    this->add_node(target);
  }

  auto source_ref = _osm_id_to_index.find(source.osm_id);
  assert(source_ref != _osm_id_to_index.end());
  auto target_ref = _osm_id_to_index.find(target.osm_id);
  assert(target_ref != _osm_id_to_index.end());

  adjacency_list[source_ref->second].emplace_back(osm_way_id,
                                                  target_ref->second,
                                                  length);
  adjacency_list[target_ref->second].emplace_back(osm_way_id,
                                                  source_ref->second,
                                                  length);
}

} // namespace posman
