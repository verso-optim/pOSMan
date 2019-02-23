#ifndef UNDIRECTED_GRAPH_H
#define UNDIRECTED_GRAPH_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <unordered_map>
#include <vector>

#include "structures/graph/node.h"

namespace posman {

struct Edge {
  Id osm_way_id;
  Index to;
  Distance length;

  Edge(Id osm_way_id, Index to, Distance length)
    : osm_way_id(osm_way_id), to(to), length(length) {
  }
};

class UndirectedGraph {
private:
  std::unordered_map<Id, Index> _osm_id_to_index;

  bool has_node(const Node& node);

  void add_node(const Node& node);

public:
  std::vector<Node> nodes;
  std::vector<std::vector<Edge>> adjacency_list;

  void add_edge(Id osm_way_id,
                const Node& source,
                const Node& target,
                Distance length);

  unsigned number_of_nodes() const;

  unsigned number_of_edges() const;

  std::vector<Distance> one_to_many(Id source, std::vector<Id> targets) const;
};

} // namespace posman

#endif
