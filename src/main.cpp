/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <iostream>
#include <unistd.h>

void display_usage() {
  std::string usage = "Copyright (C) 2019, VERSO\n";
  usage += "\n\tposman [OPTION]... -i FILE\n";
  usage += "Options:\n";
  usage += "\t-e EDGES,\t\t\t file containing OSM edges\n";
  usage += "\t-n NODES,\t\t\t file containing OSM nodes\n";
  std::cout << usage << std::endl;
  exit(0);
}

int main(int argc, char** argv) {
  // Parsing command-line arguments.
  const char* optString = "e:n:h?";
  int opt = getopt(argc, argv, optString);

  std::string edges_file;
  std::string nodes_file;

  while (opt != -1) {
    switch (opt) {
    case 'e':
      edges_file = optarg;
      break;
    case 'h':
      display_usage();
      break;
    case 'n':
      nodes_file = optarg;
      break;
    default:
      break;
    }
    opt = getopt(argc, argv, optString);
  }

  std::cout << nodes_file << " / " << edges_file << std::endl;

  return 0;
}
