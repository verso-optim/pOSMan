#ifndef MUNKRES_H
#define MUNKRES_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <unordered_map>
#include <vector>

#include "structures/typedefs.h"

namespace posman {

template <class T>
std::unordered_map<Index, Index>
minimum_weight_perfect_matching(const std::vector<std::vector<T>>& m);

template <class T>
std::unordered_map<Index, Index>
greedy_symmetric_approx_mwpm(const std::vector<std::vector<T>>& m);

} // namespace posman

#endif
