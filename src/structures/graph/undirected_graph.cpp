/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>

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
  auto source_rank = source_ref->second;
  auto target_ref = _osm_id_to_index.find(target.osm_id);
  assert(target_ref != _osm_id_to_index.end());
  auto target_rank = target_ref->second;

  assert(source.osm_id == nodes[source_rank].osm_id);
  assert(target.osm_id == nodes[target_rank].osm_id);

  adjacency_list[source_rank].emplace_back(osm_way_id, target_rank, length);
  adjacency_list[target_rank].emplace_back(osm_way_id, source_rank, length);
}

void UndirectedGraph::duplicate_edge(const Node& source, const Node& target) {
  auto source_ref = _osm_id_to_index.find(source.osm_id);
  assert(source_ref != _osm_id_to_index.end());
  auto source_rank = source_ref->second;
  auto target_ref = _osm_id_to_index.find(target.osm_id);
  assert(target_ref != _osm_id_to_index.end());
  auto target_rank = target_ref->second;

  // Find existing edge.
  auto search =
    std::find_if(adjacency_list[source_rank].begin(),
                 adjacency_list[source_rank].end(),
                 [&](const auto& e) { return e.to == target_rank; });
  assert(search != adjacency_list[source_rank].end());

  this->add_edge(search->osm_way_id, source, target, search->length);
}

unsigned UndirectedGraph::number_of_nodes() const {
  return nodes.size();
}

unsigned UndirectedGraph::number_of_edges() const {
  unsigned adjacencies =
    std::accumulate(adjacency_list.begin(),
                    adjacency_list.end(),
                    0,
                    [&](auto sum, const auto& l) { return sum + l.size(); });
  assert(adjacencies % 2 == 0);

  return adjacencies / 2;
}

unsigned UndirectedGraph::degree(Id osm_id) const {
  auto search = _osm_id_to_index.find(osm_id);
  assert(search != _osm_id_to_index.end());

  return adjacency_list[search->second].size();
}

Node UndirectedGraph::node_from_id(Id osm_id) const {
  auto search = _osm_id_to_index.find(osm_id);
  assert(search != _osm_id_to_index.end());

  return nodes[search->second];
}

std::vector<Id> UndirectedGraph::neighbours_ids(Id osm_id) const {
  auto search = _osm_id_to_index.find(osm_id);
  assert(search != _osm_id_to_index.end());

  std::vector<Id> ids;
  std::transform(adjacency_list[search->second].begin(),
                 adjacency_list[search->second].end(),
                 std::back_inserter(ids),
                 [&](const auto& e) { return nodes[e.to].osm_id; });
  return ids;
}

struct IndexWithDistance {
  Index index;
  Distance length;
  IndexWithDistance(Index index, Distance length)
    : index(index), length(length){};
  bool operator<(IndexWithDistance rhs) const {
    return this->length > rhs.length;
  }
};

std::vector<Distance>
UndirectedGraph::one_to_many(Id source,
                             const std::vector<Id>& targets,
                             std::unordered_map<Id, Id>& parents_map) const {
  std::vector<Distance> target_lengths(targets.size());

  auto source_ref = _osm_id_to_index.find(source);
  assert(source_ref != _osm_id_to_index.end());
  auto source_rank = source_ref->second;

  // Used to:
  //
  // 1. translate OSM id to index in nodes
  // 2. remember where to store lengths
  // 3. for stopping condition
  std::unordered_map<Index, std::size_t> node_rank_to_target_rank;

  for (std::size_t i = 0; i < targets.size(); ++i) {
    auto target_ref = _osm_id_to_index.find(targets[i]);
    assert(target_ref != _osm_id_to_index.end());

    node_rank_to_target_rank.insert(std::make_pair(target_ref->second, i));
  }

  // Simple Dijkstra to reach all targets.
  std::priority_queue<IndexWithDistance, std::vector<IndexWithDistance>>
    to_visit;
  to_visit.emplace(source_rank, 0);

  // Remember shortest lengths.
  std::vector<Distance> lengths(nodes.size(),
                                std::numeric_limits<Distance>::max());
  lengths[source_rank] = 0;

  while (!to_visit.empty()) {
    // Get node index with minimum length and pop it from the queue.
    auto current_node = to_visit.top();
    auto current_index = current_node.index;
    to_visit.pop();

    auto search_target = node_rank_to_target_rank.find(current_index);
    if (search_target != node_rank_to_target_rank.end()) {
      // We've reached one of the targets.
      target_lengths[search_target->second] = current_node.length;
      node_rank_to_target_rank.erase(search_target);
      if (node_rank_to_target_rank.empty()) {
        break;
      }
    }

    if (current_node.length > lengths[current_index]) {
      // Discard inserting a candidate already stored with cheaper
      // length.
      continue;
    }

    for (const auto& edge : adjacency_list[current_index]) {
      auto new_length = lengths[current_index] + edge.length;
      if (new_length < lengths[edge.to]) {
        // Update known length to node.
        lengths[edge.to] = new_length;
        to_visit.emplace(edge.to, new_length);
        // Update parent for destination node.
        auto search = parents_map.find(nodes[edge.to].osm_id);
        if (search != parents_map.end()) {
          search->second = nodes[current_index].osm_id;
        } else {
          parents_map.insert(
            std::make_pair(nodes[edge.to].osm_id, nodes[current_index].osm_id));
        }
      }
    }
  }

  if (!node_rank_to_target_rank.empty()) {
    std::cout << "Unreachable nodes: ";
    for (const auto& m : node_rank_to_target_rank) {
      std::cout << nodes[m.first].osm_id << " ; ";
    }
    std::cout << std::endl;
  }
  assert(node_rank_to_target_rank.empty());

  return target_lengths;
}

} // namespace posman
