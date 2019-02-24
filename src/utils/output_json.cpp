/*

This file is part of pOSMan.

Copyright (c) 2019, VERSO.
All rights reserved (see LICENSE).

*/

#include <algorithm>
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
      json_properties.AddMember("stroke", rapidjson::Value(), allocator);
      json_properties["stroke"].SetString(
        (current_way == META_WAY_ID)
          ? "#FF0000"
          : (current_way == META_WAY_ID - 1) ? "#0000FF" : "#555555");
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

void write_output(const UndirectedGraph& graph,
                  const std::vector<Index>& path,
                  const std::vector<Id>& path_way_ids,
                  const std::string& output_file) {
  rapidjson::Document json_output;
  json_output.SetObject();
  rapidjson::Document::AllocatorType& allocator = json_output.GetAllocator();

  rapidjson::Value json_legs(rapidjson::kArrayType);

  assert(path[0] == path[path.size() - 1]);
  assert(path_way_ids.size() + 1 == path.size());

  Distance total_length = 0;

  for (std::size_t i = 0; i < path_way_ids.size(); ++i) {
    auto current_rank = path[i];
    auto next_rank = path[i + 1];
    auto current_way_id = path_way_ids[i];

    auto edge = std::find_if(graph.adjacency_list[current_rank].begin(),
                             graph.adjacency_list[current_rank].end(),
                             [current_way_id, next_rank](const auto& e) {
                               return (e.osm_way_id == current_way_id) and
                                      (e.to == next_rank);
                             });

    assert(edge != graph.adjacency_list[current_rank].end());

    rapidjson::Value json_leg(rapidjson::kObjectType);
    json_leg.AddMember("length", edge->length, allocator);
    json_leg.AddMember("way_id", edge->osm_way_id, allocator);

    rapidjson::Value json_nodes(rapidjson::kArrayType);

    rapidjson::Value json_from_node(rapidjson::kObjectType);
    rapidjson::Value json_from_coords(rapidjson::kArrayType);
    json_from_coords.PushBack(graph.nodes[current_rank].lon, allocator);
    json_from_coords.PushBack(graph.nodes[current_rank].lat, allocator);
    json_from_node.AddMember("id", graph.nodes[current_rank].osm_id, allocator);
    json_from_node.AddMember("coordinates", json_from_coords, allocator);
    json_nodes.PushBack(json_from_node, allocator);

    rapidjson::Value json_to_node(rapidjson::kObjectType);
    rapidjson::Value json_to_coords(rapidjson::kArrayType);
    json_to_coords.PushBack(graph.nodes[next_rank].lon, allocator);
    json_to_coords.PushBack(graph.nodes[next_rank].lat, allocator);
    json_to_node.AddMember("id", graph.nodes[next_rank].osm_id, allocator);
    json_to_node.AddMember("coordinates", json_to_coords, allocator);
    json_nodes.PushBack(json_to_node, allocator);

    json_leg.AddMember("nodes", json_nodes, allocator);
    json_legs.PushBack(json_leg, allocator);

    total_length += edge->length;
  }

  json_output.AddMember("legs", json_legs, allocator);
  json_output.AddMember("length", total_length, allocator);

  write_to_json(json_output, output_file);
}

} // namespace io
} // namespace posman
