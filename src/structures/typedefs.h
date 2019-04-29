#ifndef POSMAN_TYPEDEFS_H
#define POSMAN_TYPEDEFS_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <limits>
#include <vector>

namespace posman {

using Id = uint64_t;
using Index = uint32_t;
using Distance = uint32_t;
using Coordinate = double;

using GeometryList = std::unordered_map<
  Id,
  std::unordered_map<Id, std::vector<std::array<Coordinate, 2>>>>;

constexpr Id META_WAY_ID = std::numeric_limits<Id>::max();

} // namespace posman

#endif
