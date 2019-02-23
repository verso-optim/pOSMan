/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <fstream>
#include <iostream>
#include <unordered_set>

#include "../include/rapidjson/stringbuffer.h"
#include "../include/rapidjson/writer.h"

#include "utils/output_json.h"

namespace posman {
namespace io {

inline void write_to_json(const rapidjson::Document& json_output,
                          const std::string& output_file) {
  // Rapidjson writing process.
  rapidjson::StringBuffer s;
  rapidjson::Writer<rapidjson::StringBuffer> r_writer(s);
  json_output.Accept(r_writer);

  // Write to relevant output.
  if (output_file.empty()) {
    // Log to standard output.
    std::cout << s.GetString() << std::endl;
  } else {
    // Log to file.
    std::ofstream out_stream(output_file, std::ofstream::out);
    out_stream << s.GetString();
    out_stream.close();
  }
}

void log_graph_as_geojson(const UndirectedGraph& graph,
                          const std::string& output_file) {
  rapidjson::Document json_output;
  json_output.SetObject();
  rapidjson::Document::AllocatorType& allocator = json_output.GetAllocator();

  json_output.AddMember("type", rapidjson::Value(), allocator);
  json_output["type"].SetString("FeatureCollection");

  rapidjson::Value json_features(rapidjson::kArrayType);

  for (std::size_t i = 0; i < graph.number_of_nodes(); ++i) {
    for (const auto& edge : graph.adjacency_list[i]) {
      auto current_way = edge.osm_way_id;

      rapidjson::Value json_feature(rapidjson::kObjectType);
      json_feature.AddMember("type", rapidjson::Value(), allocator);
      json_feature["type"].SetString("Feature");

      rapidjson::Value json_properties(rapidjson::kObjectType);
      json_properties.AddMember("way", current_way, allocator);
      json_properties.AddMember("source", graph.nodes[i].osm_id, allocator);
      json_properties.AddMember("target",
                                graph.nodes[edge.to].osm_id,
                                allocator);
      json_feature.AddMember("properties", json_properties, allocator);

      rapidjson::Value json_geometry(rapidjson::kObjectType);
      json_geometry.AddMember("type", rapidjson::Value(), allocator);
      json_geometry["type"].SetString("LineString");

      rapidjson::Value json_coordinates(rapidjson::kArrayType);
      rapidjson::Value json_source(rapidjson::kArrayType);
      json_source.PushBack(graph.nodes[i].lon, allocator);
      json_source.PushBack(graph.nodes[i].lat, allocator);
      json_coordinates.PushBack(json_source, allocator);
      rapidjson::Value json_target(rapidjson::kArrayType);
      json_target.PushBack(graph.nodes[edge.to].lon, allocator);
      json_target.PushBack(graph.nodes[edge.to].lat, allocator);
      json_coordinates.PushBack(json_target, allocator);
      json_geometry.AddMember("coordinates", json_coordinates, allocator);

      json_feature.AddMember("geometry", json_geometry, allocator);

      json_features.PushBack(json_feature, allocator);
    }
  }

  json_output.AddMember("features", json_features, allocator);

  write_to_json(json_output, output_file);
}

} // namespace io
} // namespace posman
