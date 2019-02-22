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

  Node(Id osm_id, Coordinate lon, Coordinate lat)
    : osm_id(osm_id), lon(lon), lat(lat) {
  }

  Node(const std::string& flat) {
    // id
    std::string::size_type e = flat.find(',');
    osm_id = strtoul(flat.substr(0, e).c_str(), nullptr, 10);
    // lon
    std::string::size_type s = e + 1;
    e = flat.find(',', s);
    lon = std::stod(flat.substr(s, e - s));
    // lat
    s = e + 1;
    e = flat.size();
    lat = std::stod(flat.substr(s, e - s));
  }
};

} // namespace posman

#endif
