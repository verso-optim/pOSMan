#ifndef MATCHING_H
#define MATCHING_H

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
compute_matching(const std::vector<std::vector<T>>& m);

template <class T>
std::unordered_map<Index, Index>
blossom_matching(const std::vector<std::vector<T>>& m);

} // namespace posman

#endif
