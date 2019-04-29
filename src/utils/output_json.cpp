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

inline void
get_detailed_geometry(const GeometryList& geometries,
                      const Node& source_node,
                      const Node& target_node,
                      rapidjson::Value& json_geometry,
                      rapidjson::Document::AllocatorType& allocator) {
  auto source_id = source_node.osm_id;
  auto target_id = target_node.osm_id;
  bool reverse = false;

  auto search = geometries.find(source_id);
  auto geom = ((geometries.begin())->second).begin();

  if (search == geometries.end()) {
    reverse = true;
  } else {
    geom = (search->second).find(target_id);
    if (geom == (search->second).end()) {
      reverse = true;
    }
  }

  if (reverse) {
    search = geometries.find(target_id);
    assert(search != geometries.end());
    geom = (search->second).find(source_id);
    assert(geom != (search->second).end());
  }

  auto size = (geom->second).size();
  for (unsigned i = 0; i < size; ++i) {
    auto& coord = (reverse) ? (geom->second)[size - 1 - i] : (geom->second)[i];
    rapidjson::Value json_coord(rapidjson::kArrayType);
    json_coord.PushBack(coord[0], allocator);
    json_coord.PushBack(coord[1], allocator);
    json_geometry.PushBack(json_coord, allocator);
  }
}

void log_graph_as_geojson(const UndirectedGraph& graph,
                          const std::string& output_file,
                          const GeometryList& geometries) {
  rapidjson::Document json_output;
  json_output.SetObject();
  rapidjson::Document::AllocatorType& allocator = json_output.GetAllocator();

  json_output.AddMember("type", rapidjson::Value(), allocator);
  json_output["type"].SetString("FeatureCollection");

  rapidjson::Value json_features(rapidjson::kArrayType);

  for (std::size_t i = 0; i < graph.number_of_nodes(); ++i) {
    const auto& node = graph.nodes[i];
    for (const auto& edge : graph.adjacency_list[i]) {
      auto current_way = edge.osm_way_id;

      rapidjson::Value json_feature(rapidjson::kObjectType);
      json_feature.AddMember("type", rapidjson::Value(), allocator);
      json_feature["type"].SetString("Feature");

      rapidjson::Value json_properties(rapidjson::kObjectType);
      json_properties.AddMember("way", current_way, allocator);
      json_properties.AddMember("source", node.osm_id, allocator);
      json_properties.AddMember("target",
                                graph.nodes[edge.to].osm_id,
                                allocator);
      json_properties.AddMember("stroke", rapidjson::Value(), allocator);
      json_properties["stroke"].SetString(
        (current_way == META_WAY_ID)
          ? "#FF0000"
          : (current_way == META_WAY_ID - 1) ? "#0000FF" : "#555555");
      json_feature.AddMember("properties", json_properties, allocator);

      rapidjson::Value geojson_geometry(rapidjson::kObjectType);
      geojson_geometry.AddMember("type", rapidjson::Value(), allocator);
      geojson_geometry["type"].SetString("LineString");

      rapidjson::Value json_geometry(rapidjson::kArrayType);
      if (current_way == META_WAY_ID) {
        rapidjson::Value json_source(rapidjson::kArrayType);
        json_source.PushBack(node.lon, allocator);
        json_source.PushBack(node.lat, allocator);
        json_geometry.PushBack(json_source, allocator);
        rapidjson::Value json_target(rapidjson::kArrayType);
        json_target.PushBack(graph.nodes[edge.to].lon, allocator);
        json_target.PushBack(graph.nodes[edge.to].lat, allocator);
        json_geometry.PushBack(json_target, allocator);
      } else {
        // Retrieve detailed geometry.
        get_detailed_geometry(geometries,
                              node,
                              graph.nodes[edge.to],
                              json_geometry,
                              allocator);
      }

      geojson_geometry.AddMember("coordinates", json_geometry, allocator);

      json_feature.AddMember("geometry", geojson_geometry, allocator);

      json_features.PushBack(json_feature, allocator);
    }
  }

  json_output.AddMember("features", json_features, allocator);

  write_to_json(json_output, output_file);
}

void write_output(
  const UndirectedGraph& global_graph,
  const UndirectedGraph& target_graph,
  const std::vector<Index>& path,
  const std::vector<Id>& path_way_ids,
  const std::unordered_map<Id, std::unordered_map<Id, Id>>& source_to_parents,
  const GeometryList& geometries,
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

    auto edge = std::find_if(target_graph.adjacency_list[current_rank].begin(),
                             target_graph.adjacency_list[current_rank].end(),
                             [current_way_id, next_rank](const auto& e) {
                               return (e.osm_way_id == current_way_id) and
                                      (e.to == next_rank);
                             });

    assert(edge != target_graph.adjacency_list[current_rank].end());

    rapidjson::Value json_leg(rapidjson::kObjectType);
    json_leg.AddMember("length", edge->length, allocator);

    if (edge->osm_way_id < META_WAY_ID - 1) {
      // Regular edge from an OSM way.
      json_leg.AddMember("way_id", edge->osm_way_id, allocator);

      rapidjson::Value json_nodes(rapidjson::kArrayType);

      rapidjson::Value json_from_node(rapidjson::kObjectType);
      rapidjson::Value json_from_coords(rapidjson::kArrayType);
      json_from_coords.PushBack(target_graph.nodes[current_rank].lon,
                                allocator);
      json_from_coords.PushBack(target_graph.nodes[current_rank].lat,
                                allocator);
      json_from_node.AddMember("id",
                               target_graph.nodes[current_rank].osm_id,
                               allocator);
      json_from_node.AddMember("coordinates", json_from_coords, allocator);
      json_nodes.PushBack(json_from_node, allocator);

      rapidjson::Value json_to_node(rapidjson::kObjectType);
      rapidjson::Value json_to_coords(rapidjson::kArrayType);
      json_to_coords.PushBack(target_graph.nodes[next_rank].lon, allocator);
      json_to_coords.PushBack(target_graph.nodes[next_rank].lat, allocator);
      json_to_node.AddMember("id",
                             target_graph.nodes[next_rank].osm_id,
                             allocator);
      json_to_node.AddMember("coordinates", json_to_coords, allocator);
      json_nodes.PushBack(json_to_node, allocator);

      json_leg.AddMember("nodes", json_nodes, allocator);

      rapidjson::Value json_geometry(rapidjson::kArrayType);
      get_detailed_geometry(geometries,
                            target_graph.nodes[current_rank],
                            target_graph.nodes[next_rank],
                            json_geometry,
                            allocator);

      json_leg.AddMember("geometry", json_geometry, allocator);
      json_legs.PushBack(json_leg, allocator);
    } else {
      // Meta edge added after one-to-many search, we want to retrieve
      // the whole path.
      json_leg.AddMember("way_id", 0, allocator);

      auto current_id = target_graph.nodes[current_rank].osm_id;
      auto search_map = source_to_parents.find(current_id);
      assert(search_map != source_to_parents.end());

      auto& parents_map = search_map->second;

      std::vector<Id> reversed_ids;
      Id dijkstra_id = target_graph.nodes[next_rank].osm_id;
      while (dijkstra_id != current_id) {
        reversed_ids.push_back(dijkstra_id);
        auto search = parents_map.find(dijkstra_id);
        assert(search != parents_map.end());
        dijkstra_id = search->second;
      }
      reversed_ids.push_back(current_id);

      rapidjson::Value json_nodes(rapidjson::kArrayType);
      rapidjson::Value json_geometry(rapidjson::kArrayType);
      for (auto iter = reversed_ids.crbegin(); iter != reversed_ids.crend();
           ++iter) {
        const auto& current_node = global_graph.node_from_id(*iter);
        assert(current_node.osm_id == *iter);

        auto next = iter + 1;
        if (next != reversed_ids.crend()) {
          const auto& next_node = global_graph.node_from_id(*(next));
          if (json_geometry.Size() > 0) {
            json_geometry.PopBack();
          }
          get_detailed_geometry(geometries,
                                current_node,
                                next_node,
                                json_geometry,
                                allocator);
        }

        rapidjson::Value json_node(rapidjson::kObjectType);
        rapidjson::Value json_coords(rapidjson::kArrayType);
        json_coords.PushBack(current_node.lon, allocator);
        json_coords.PushBack(current_node.lat, allocator);
        json_node.AddMember("id", current_node.osm_id, allocator);
        json_node.AddMember("coordinates", json_coords, allocator);
        json_nodes.PushBack(json_node, allocator);
      }

      json_leg.AddMember("nodes", json_nodes, allocator);
      json_leg.AddMember("geometry", json_geometry, allocator);
      json_legs.PushBack(json_leg, allocator);
    }

    total_length += edge->length;
  }

  json_output.AddMember("legs", json_legs, allocator);
  json_output.AddMember("length", total_length, allocator);

  write_to_json(json_output, output_file);
}

} // namespace io
} // namespace posman
