/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unistd.h>

#include "algorithms/munkres.h"
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

  using namespace posman;

  // 1. Parsing the whole graph from OSM file.
  std::cout << "[info] Loading global graph." << std::endl;
  auto global_graph = io::parse_graph(nodes_file, edges_file);

  std::cout << "  Global graph has " << global_graph.number_of_nodes()
            << " nodes and " << global_graph.number_of_edges() << " edges."
            << std::endl;

  // 2. Creating the target graph based on provided OSM way ids.
  std::cout << "[info] Loading target graph." << std::endl;
  auto target_graph = io::parse_ways(ways_file, global_graph);

  unsigned target_nodes = target_graph.number_of_nodes();

  std::cout << "  Target graph has " << target_nodes << " nodes and "
            << target_graph.number_of_edges() << " edges." << std::endl;

  // 3. Find out all nodes of odd degree in target rank.
  std::vector<Index> odd_degree_node_ranks;
  for (std::size_t i = 0; i < target_nodes; ++i) {
    if (target_graph.adjacency_list[i].size() % 2 != 0) {
      odd_degree_node_ranks.push_back(i);
    }
  }

  auto odd_nodes = odd_degree_node_ranks.size();
  std::cout << "  Target graph has " << odd_nodes << " nodes with odd degree ("
            << static_cast<int>(100 * static_cast<float>(odd_nodes) /
                                target_nodes)
            << "%)" << std::endl;

  // 4. Compute a many-to-many matrix of distances (in the global
  // graph) between nodes of odd degree in the target graph.
  std::vector<Id> node_ids;
  std::transform(odd_degree_node_ranks.begin(),
                 odd_degree_node_ranks.end(),
                 std::back_inserter(node_ids),
                 [&](auto rank) { return target_graph.nodes[rank].osm_id; });

  auto N = node_ids.size();
  std::vector<std::vector<Distance>> lengths_matrix(N,
                                                    std::vector<Distance>(N));

  for (std::size_t i = 0; i < N; ++i) {
    auto line = global_graph.one_to_many(node_ids[0], node_ids);
    for (std::size_t j = 0; j < line.size(); ++j) {
      lengths_matrix[i][i + j] = line[j];
      lengths_matrix[i + j][i] = line[j];
    }
    lengths_matrix[i][i] = 3 * (std::numeric_limits<Distance>::max() / 4);
    node_ids.erase(node_ids.begin());
  }

  // 5. Compute minimum weight perfect matching for odd nodes.
  auto mwpm = minimum_weight_perfect_matching(lengths_matrix);

  // Storing those edges from mwpm that are coherent regarding
  // symmetry (y -> x whenever x -> y). Remembering the rest of them
  // for further use.
  std::vector<Index> wrong_node_ranks;
  unsigned total_ok = 0;

  for (const auto& match : mwpm) {
    assert(match.first != match.second);

    auto first_node_rank = odd_degree_node_ranks[match.first];
    auto second_node_rank = odd_degree_node_ranks[match.second];

    if (mwpm.at(match.second) == match.first) {
      ++total_ok;
      target_graph.add_edge(META_WAY_ID,
                            target_graph.nodes[first_node_rank],
                            target_graph.nodes[second_node_rank],
                            lengths_matrix[match.first][match.second]);
    } else {
      wrong_node_ranks.push_back(match.first);
    }
  }

  // Log the target graph for easy visualisation.
  if (!geojson_target.empty()) {
    io::log_graph_as_geojson(target_graph, geojson_target);
  }

  return 0;
}
