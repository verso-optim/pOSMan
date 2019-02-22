#ifndef NODE_H
#define NODE_H

/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include "structures/typedefs.h"

namespace posman {

struct Node {
  Id osm_id;
  Coordinate lon;
  Coordinate lat;
};

} // namespace posman

#endif
