/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_set>

#include "algorithms/matching.h"
#include "algorithms/munkres.h"

namespace posman {

template <class T>
std::unordered_map<Index, Index>
compute_matching(const std::vector<std::vector<T>>& m) {
  std::unordered_map<Index, Index> matching;

  auto mwpm = minimum_weight_perfect_matching(m);

  // Storing those edges from mwpm that are coherent regarding
  // symmetry (y -> x whenever x -> y). Remembering the rest of them
  // by their rank in odd_degree_node_ranks for further use.
  std::vector<Index> wrong_node_ranks;

  std::unordered_set<Index> already_added;
  for (const auto& match : mwpm) {
    assert(match.first != match.second);

    if (mwpm[match.second] == match.first) {
      // Real symmetric edge.
      if (already_added.find(match.first) != already_added.end()) {
        continue;
      }
      already_added.insert(match.second);

      matching.emplace(match.first, match.second);
    } else {
      wrong_node_ranks.push_back(match.first);
    }
  }

  if (!wrong_node_ranks.empty()) {
    // Spot cycles.
    std::unordered_set<Index> seen;
    std::vector<std::vector<Index>> cycles;
    std::vector<std::vector<T>> cycle_costs;

    for (const auto wn : wrong_node_ranks) {
      if (seen.find(wn) != seen.end()) {
        continue;
      }
      seen.insert(wn);
      cycles.push_back({wn});
      auto& current_cycle = cycles.back();
      cycle_costs.push_back({m[wn][mwpm[wn]]});
      auto& current_costs = cycle_costs.back();

      auto current = wn;
      do {
        auto next = mwpm[current];
        seen.insert(next);

        current_cycle.push_back(next);
        current_costs.push_back(m[next][mwpm[next]]);

        current = next;
      } while (mwpm[current] != wn);
    }

    for (std::size_t i = 0; i < cycles.size(); ++i) {
      auto& cycle = cycles[i];
      if (cycle.size() % 2 == 0) {
        std::cout << "[Info] Force even cycle: ";
        for (std::size_t k = 0; k < cycle.size(); ++k) {
          std::cout << cycle[k] << "\t";
          auto search = std::find(wrong_node_ranks.begin(),
                                  wrong_node_ranks.end(),
                                  cycle[k]);
          assert(search != wrong_node_ranks.end());
          wrong_node_ranks.erase(search);

          if (k % 2 == 0) {
            auto first_rank = cycle[k];
            auto second_rank = cycle[(k + 1) % cycle.size()];
            matching.emplace(first_rank, second_rank);
          }
        }
        std::cout << std::endl;
      } else {
        auto& costs = cycle_costs[i];

        auto min_search = std::min_element(costs.begin(), costs.end());
        auto offset = std::distance(costs.begin(), min_search);
        auto first_rank = cycle[offset];
        auto second_rank = cycle[(offset + 1) % cycle.size()];

        matching.emplace(first_rank, second_rank);

        std::cout << "[Info] Pick cheapest edge in odd cycle: " << first_rank
                  << " -> " << second_rank << std::endl;

        auto search = std::find(wrong_node_ranks.begin(),
                                wrong_node_ranks.end(),
                                first_rank);
        assert(search != wrong_node_ranks.end());
        wrong_node_ranks.erase(search);
        search = std::find(wrong_node_ranks.begin(),
                           wrong_node_ranks.end(),
                           second_rank);
        assert(search != wrong_node_ranks.end());
        wrong_node_ranks.erase(search);
      }
    }

    // Run a greedy algorithm to get a symmetric matching for wrong
    // nodes.
    std::vector<std::vector<T>> sub_matrix(wrong_node_ranks.size(),
                                           std::vector<T>(
                                             wrong_node_ranks.size()));
    for (std::size_t i = 0; i < wrong_node_ranks.size(); ++i) {
      for (std::size_t j = 0; j < wrong_node_ranks.size(); ++j) {
        sub_matrix[i][j] = m[wrong_node_ranks[i]][wrong_node_ranks[j]];
      }
    }

    auto remaining_greedy_mwpm = greedy_symmetric_approx_mwpm(sub_matrix);

    for (const auto& match : remaining_greedy_mwpm) {
      assert(matching.find(wrong_node_ranks[match.second]) == matching.end());
      matching.emplace(wrong_node_ranks[match.first],
                       wrong_node_ranks[match.second]);
    }
  }

  return matching;
}

template std::unordered_map<Index, Index>
compute_matching(const std::vector<std::vector<Distance>>& m);

} // namespace posman
