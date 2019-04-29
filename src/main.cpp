/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unistd.h>
#include <unordered_set>

#include "algorithms/matching.h"
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
  usage += "\t-o OUTPUT,\t\t\t output file name\n";
  usage += "\t-s START,\t\t\t OSM node id for starting point\n";
  usage += "\t-w WAYS,\t\t\t file containing target ways\n";
  std::cout << usage << std::endl;

  exit(0);
}

int main(int argc, char** argv) {
  // Parsing command-line arguments.
  const char* optString = "e:g:h?n:o:s:w:";
  int opt = getopt(argc, argv, optString);

  std::string edges_file;
  std::string geojson_target;
  std::string nodes_file;
  std::string output_file;
  posman::Id start_id;
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
    case 'o':
      output_file = optarg;
      break;
    case 's':
      start_id = std::stoul(optarg);
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
  GeometryList geometries;
  auto global_graph = io::parse_graph(nodes_file, edges_file, geometries);

  std::cout << "  Global graph has " << global_graph.number_of_nodes()
            << " nodes and " << global_graph.number_of_edges() << " edges."
            << std::endl;

  // 2. Creating the target graph based on provided OSM way ids.
  std::cout << "[info] Loading target graph." << std::endl;
  auto target_graph = io::parse_ways(ways_file, global_graph);

  unsigned target_nodes = target_graph.number_of_nodes();

  std::cout << "  Target graph has " << target_nodes << " nodes and "
            << target_graph.number_of_edges() << " edges." << std::endl;

  // 2. bis Force doubling edges for dead-ends because this is what we
  // would like to do for those in step 3. anyway (and it avoids
  // potential non-symmetric matching in Munkres).
  std::vector<std::vector<Node>> addition_chains;

  for (std::size_t i = 0; i < target_nodes; ++i) {
    Node current = target_graph.nodes[i];
    if (global_graph.degree(current.osm_id) == 1) {
      auto neighbours_ids = target_graph.neighbours_ids(current.osm_id);
      assert(neighbours_ids.size() == 1);

      Node next = target_graph.node_from_id(neighbours_ids[0]);

      addition_chains.push_back({current, next});
      auto& addition_chain = addition_chains.back();

      auto previous_id = current.osm_id;
      current = next;
      while (target_graph.degree(current.osm_id) == 2) {
        neighbours_ids = target_graph.neighbours_ids(current.osm_id);
        assert(neighbours_ids.size() == 2);

        Id next_id = (neighbours_ids[0] == previous_id) ? neighbours_ids[1]
                                                        : neighbours_ids[0];
        assert(next_id != previous_id);

        next = target_graph.node_from_id(next_id);
        addition_chain.push_back(next);

        previous_id = current.osm_id;
        current = next;
      }
    }
  }

  for (const auto& addition_chain : addition_chains) {
    for (unsigned i = 0; i < addition_chain.size() - 1; ++i) {
      target_graph.duplicate_edge(addition_chain[i], addition_chain[i + 1]);
    }
  }

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
  std::vector<Id> odd_nodes_ids;
  std::transform(odd_degree_node_ranks.begin(),
                 odd_degree_node_ranks.end(),
                 std::back_inserter(odd_nodes_ids),
                 [&](const auto rank) {
                   return target_graph.nodes[rank].osm_id;
                 });

  auto N = odd_nodes_ids.size();
  std::vector<std::vector<Distance>> lengths_matrix(N);

  //  Remember parent nodes for each source in one-to-many search.
  std::unordered_map<Id, std::unordered_map<Id, Id>> source_to_parents;

  for (std::size_t i = 0; i < N; ++i) {
    auto source_map = source_to_parents.insert(
      std::make_pair(target_graph.nodes[odd_degree_node_ranks[i]].osm_id,
                     std::unordered_map<Id, Id>()));

    assert(source_map.second);

    lengths_matrix[i] = global_graph.one_to_many(odd_nodes_ids[i],
                                                 odd_nodes_ids,
                                                 source_map.first->second);
    lengths_matrix[i][i] = 3 * (std::numeric_limits<Distance>::max() / 4);
  }

  // 5. Compute perfect matching for odd nodes.
  auto matching = compute_matching(lengths_matrix);

  assert(2 * matching.size() == odd_degree_node_ranks.size());

  for (const auto& match : matching) {
    auto first_node_rank = odd_degree_node_ranks[match.first];
    auto second_node_rank = odd_degree_node_ranks[match.second];

    target_graph.add_edge(META_WAY_ID,
                          target_graph.nodes[first_node_rank],
                          target_graph.nodes[second_node_rank],
                          lengths_matrix[match.first][match.second]);
  }

  // Target graph is now eulerian, log for easy visualisation.
  assert(std::count_if(target_graph.adjacency_list.begin(),
                       target_graph.adjacency_list.end(),
                       [](const auto& l) { return l.size() % 2 != 0; }) == 0);

  if (!geojson_target.empty()) {
    io::log_graph_as_geojson(target_graph, geojson_target, geometries);
  }

  // 6. Use Hierholzer's algorithm to derive an eulerian path.
  auto eulerian_graph = target_graph;
  auto start =
    std::find_if(eulerian_graph.nodes.begin(),
                 eulerian_graph.nodes.end(),
                 [&](const auto& n) { return n.osm_id == start_id; });
  assert(start != eulerian_graph.nodes.end());

  std::vector<Index> eulerian_path;
  std::vector<Id> eulerian_path_way_ids;
  eulerian_path.push_back(std::distance(eulerian_graph.nodes.begin(), start));

  bool complete_tour = false;
  while (!complete_tour) {
    // Finding first element of eulerian_path that still has an
    // adjacent edge (if any).
    auto new_tour_start =
      std::find_if(eulerian_path.begin(),
                   eulerian_path.end(),
                   [&](const auto rank) {
                     return !eulerian_graph.adjacency_list[rank].empty();
                   });
    complete_tour = (new_tour_start == eulerian_path.end());

    if (!complete_tour) {
      // Add new tour to initial eulerian path and check again.
      std::vector<Index> new_tour;
      std::vector<Id> new_ways;
      Index initial_node_rank = *new_tour_start;
      Index current_node_rank = initial_node_rank;
      Index next_node_rank;
      do {
        new_tour.push_back(current_node_rank);
        // Find next node from any adjacent edge and remove used edge
        // in both ways in adjacency lists.
        auto chosen_edge =
          eulerian_graph.adjacency_list[current_node_rank].begin();
        next_node_rank = chosen_edge->to;
        new_ways.push_back(chosen_edge->osm_way_id);

        auto reverse_edge =
          std::find_if(eulerian_graph.adjacency_list[next_node_rank].begin(),
                       eulerian_graph.adjacency_list[next_node_rank].end(),
                       [&](const auto& e) {
                         return (e.osm_way_id == chosen_edge->osm_way_id) and
                                (e.to == current_node_rank);
                       });

        assert(reverse_edge !=
               eulerian_graph.adjacency_list[next_node_rank].end());
        eulerian_graph.adjacency_list[current_node_rank].erase(chosen_edge);
        eulerian_graph.adjacency_list[next_node_rank].erase(reverse_edge);

        current_node_rank = next_node_rank;
      } while (current_node_rank != initial_node_rank);

      // Adding new tour to existing eulerian path.
      auto new_ways_start =
        eulerian_path_way_ids.begin() +
        std::distance(eulerian_path.begin(), new_tour_start);

      eulerian_path.insert(new_tour_start, new_tour.begin(), new_tour.end());
      eulerian_path_way_ids.insert(new_ways_start,
                                   new_ways.begin(),
                                   new_ways.end());
    }
  }

  if (!output_file.empty()) {
    io::write_output(global_graph,
                     target_graph,
                     eulerian_path,
                     eulerian_path_way_ids,
                     source_to_parents,
                     geometries,
                     output_file);
  }

  return 0;
}
