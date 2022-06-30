//
// Created by sunsx on 28/05/21.
//

#include <queue>
#include <computesetintersection.h>
#include "order_manager.h"
#include "graphoperations.h"

void OrderManager::initialize(const Graph *query_graph) {}

void OrderManager::release() {
    automorphisms_.clear();
    orders_.clear();
    edge_orders_mapping_.clear();
    label_edge_mapping_.clear();
    label_automorphism_mapping_.clear();
    automorphism_edges_.clear();
}

OrdersPerEdge *OrderManager::get_orders(Edge edge) {
    auto it = edge_orders_mapping_.find(edge);
    if (it != edge_orders_mapping_.end()) {
        return &orders_[it->second];
    }

    return nullptr;
}

std::vector<uint32_t> *OrderManager::get_mapped_automorphism(LabelTriple label_triple) {
    auto it = label_automorphism_mapping_.find(label_triple);
    if (it != label_automorphism_mapping_.end()) {
        return &it->second;
    }

    return nullptr;
}

std::vector<Edge> *OrderManager::get_mapped_edges(LabelTriple label_triple) {
    auto it = label_edge_mapping_.find(label_triple);
    if (it != label_edge_mapping_.end()) {
        return &it->second;
    }
    return nullptr;
}

std::tuple<Edge, OrdersPerEdge*, uint32_t> OrderManager::get_automorphism_meta(uint32_t id) {
    return std::make_tuple(automorphism_edges_[id][0], &orders_[id], (uint32_t)(automorphism_edges_[id].size()));
}

void OrderManager::detect_automorphism_edges(const Graph *query_graph) {
    // Divide the vertex into disjoint sets based on automorphisms.
    GraphOperations::compute_automorphism(query_graph, automorphisms_);

    uint32_t n = query_graph->getVerticesCount();
    spp::sparse_hash_set<Edge> selected;

    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_nbr_count;
        auto u_nbr = query_graph->getVertexNeighbors(u, u_nbr_count);

        for (uint32_t i = 0; i < u_nbr_count; ++i) {
            uint32_t uu = u_nbr[i];
            Edge e = {u, uu};
            if (!selected.contains(e)) {
                selected.insert(e);
                automorphism_edges_.push_back({e});
                for (auto &embedding: automorphisms_) {
                    Edge mapped_e = {embedding[u], embedding[uu]};
                    if (!selected.contains(mapped_e)) {
                        selected.insert(mapped_e);
                        automorphism_edges_.back().push_back(mapped_e);
                    }
                }
            }
        }
    }
}

void OrderManager::create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    uint32_t n = query_graph->getVerticesCount();
    auto& indexing_order = order.indexing_order_;
    auto& indexing_order_bn = order.indexing_order_bn_;
    auto& indexing_order_bn_offset = order.indexing_order_bn_offset_;

    indexing_order = {edge.first, edge.second};
    indexing_order_bn = {edge.first};
    indexing_order_bn_offset = {0, 0, 1};

    std::vector<bool> visited(n);
    visited[edge.first] = true;
    visited[edge.second] = true;

    std::vector<uint32_t> adjacent_to_both;
    std::vector<uint32_t> adjacent_to_one;
    std::vector<uint32_t> adjacent_to_none;

    for (uint32_t u = 0; u < n; ++u) {
        if (u == edge.first || u == edge.second)
            continue;

        bool f1 = query_graph->checkEdgeExistence(u, edge.first);
        bool f2 = query_graph->checkEdgeExistence(u, edge.second);

        if (f1 && f2) {
            adjacent_to_both.push_back(u);
        }
        else if (!f1 && !f2) {
            adjacent_to_none.push_back(u);
        }
        else {
            adjacent_to_one.push_back(u);
        }
    }

    auto generate_function = [query_graph, &visited, &indexing_order, &indexing_order_bn, &indexing_order_bn_offset]
            (std::vector<uint32_t>& target_vertex) {
        for (uint32_t i = 0; i < target_vertex.size(); ++i) {
            std::vector<uint32_t> current_vertex_bn;
            std::vector<uint32_t> selected_vertex_bn;
            uint32_t selected_vertex;

            for (auto u : target_vertex) {
                if (!visited[u]) {
                    current_vertex_bn.clear();
                    // Get backward neighbors.
                    for (auto uu : indexing_order) {
                        if (query_graph->checkEdgeExistence(u, uu)) {
                            current_vertex_bn.push_back(uu);
                        }
                    }
                    if (current_vertex_bn.size() > selected_vertex_bn.size()) {
                        current_vertex_bn.swap(selected_vertex_bn);
                        selected_vertex = u;
                    }
                }
            }

            indexing_order.push_back(selected_vertex);
            indexing_order_bn.insert(indexing_order_bn.end(), selected_vertex_bn.begin(), selected_vertex_bn.end());
            indexing_order_bn_offset.push_back(indexing_order_bn.size());

            visited[selected_vertex] = true;
        }
    };
    order.triangle_end_index_ = adjacent_to_both.size() + 2;
    order.adjacent_end_index_ = adjacent_to_both.size() + adjacent_to_one.size() + 2;

    generate_function(adjacent_to_both);
    generate_function(adjacent_to_one);
    generate_function(adjacent_to_none);
}

void OrderManager::create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order) {
    /**
     * 1. Get all connected components.
     * 2. Generate the matching order for each connected component.
     * 3. Merge the connected component.
     */
     uint32_t n = query_graph->getVerticesCount();
     std::vector<bool> visited(n, false);
     visited[edge.first] = true;
     visited[edge.second] = true;

     std::vector<std::vector<uint32_t>> connected_components;
     std::queue<uint32_t> q;

     // 1. Generate connected components by conducting a BFS from each vertex.
     for (uint32_t u = 0; u < n; ++u) {
         if (!visited[u]) {
             std::vector<uint32_t> component;
             component.push_back(u);
             q.push(u);
             visited[u] = true;

             while (!q.empty()) {
                 uint32_t uu = q.front();
                 q.pop();
                 uint32_t uu_nbrs_count;
                 auto uu_nbrs = query_graph->getVertexNeighbors(uu, uu_nbrs_count);

                 for (uint32_t i = 0; i < uu_nbrs_count; ++i) {
                     uint32_t uuu = uu_nbrs[i];
                     if (!visited[uuu]) {
                         q.push(uuu);
                         component.push_back(uuu);
                         visited[uuu] = true;
                     }
                 }
             }

             connected_components.emplace_back(component);
         }
     }

     // 2. Generate the matching order for each connected component.
     std::vector<std::vector<uint32_t>> matching_orders;

     // map the old vertex id to new vertex id
     std::vector<uint32_t> reverse_mapping(n, 0);

     std::vector<Graph*> graphs;

     for (auto& component : connected_components) {
         if (component.size() == 1 || component.size() == 2) {
             matching_orders.push_back(component);
             graphs.push_back(nullptr);
         }
         else {
             for (int i = 0; i < component.size(); ++i) {
                reverse_mapping[component[i]] = i;
             }

             // Create vertex list (vertex_id, vertex_label_id)  and edge list (begin_vertex_id, end_vertex_id, edge_label_id)
             std::vector<std::pair<uint32_t, uint32_t>> vertex_list;
             std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edge_list;

             for (auto u: component) {
                 uint32_t u_id = reverse_mapping[u];
                 uint32_t u_label = query_graph->getVertexLabel(u);
                 vertex_list.emplace_back(u_id, u_label);

                 uint32_t u_nbrs_count;
                 auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

                 for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                     uint32_t uu = u_nbrs[i];
                     if (uu != edge.first && uu != edge.second) {
                         uint32_t uu_id = reverse_mapping[uu];

                         if (u_id < uu_id) {
                             edge_list.emplace_back(u_id, uu_id, query_graph->getEdgeLabelByVertex(u_id, uu_id));
                         }
                     }
                 }
             }

             Graph *graph = new Graph(false);
             graph->is_edge_labeled = true;
             graph->loadGraphFromMemory(vertex_list, edge_list);

             std::vector<uint32_t> matching_order;
             generate_matching_order_with_RI(graph, matching_order);
             matching_orders.push_back(matching_order);
             graphs.push_back(graph);
         }
     }

     // 3. Merge the matching orders into one based on 1) the 2-core size; and 2) the component size.
     auto& updated_matching_order = order.matching_order_;
     auto& updated_matching_order_bn = order.matching_order_bn_;
     auto& updated_matching_order_bn_offset = order.matching_order_bn_offset_;
     auto& updated_matching_order_view_mapping = order.matching_order_view_mappings_;
     auto& updated_matching_order_edge_type = order.matching_order_edge_type_;

     updated_matching_order_bn_offset = {0};

     std::vector<uint32_t> core_size;
     std::vector<uint32_t> graph_size;

     // Initialize the core size and graph size for each subgraph.
     for (uint32_t i = 0; i < graphs.size(); ++i) {
         Graph* graph = graphs[i];
         if (graph != nullptr) {
             uint32_t size = 0;
             std::vector<int> core(graph->getVerticesCount(), 0);
             GraphOperations::getKCore(graph, core.data());

             for (auto core_value: core) {
                 if (core_value >= 2) {
                     size += 1;
                 }
             }
             core_size.push_back(size);
             graph_size.push_back(graph->getVerticesCount());
         }
         else {
             core_size.push_back(0);
             graph_size.push_back(connected_components[i].size());
         }
     }

     std::fill(visited.begin(), visited.end(), false);

     for (uint32_t i = 0; i < graphs.size(); ++i) {
         // Pick the graph.
         uint32_t selected_core_value = 0;
         uint32_t selected_graph_size = 0;
         uint32_t selected_graph = 0;

         for (uint32_t j = 0; j < graphs.size(); ++j) {
             if (!visited[j]) {
                 if ((core_size[j] > selected_core_value)
                     || (core_size[j] == selected_core_value && graph_size[j] > selected_graph_size)) {
                     selected_graph = j;
                     selected_core_value = core_size[j];
                     selected_graph_size = graph_size[j];
                 }
             }
         }

         visited[selected_graph] = true;

         // Update matching order and backward neighbors.
         if (graph_size[selected_graph] == 1) {
             uint32_t u = connected_components[selected_graph][0];
             updated_matching_order.push_back(u);
             if (query_graph->getVertexDegree(u) == 1) {
                 if (query_graph->checkEdgeExistence(u, edge.first)) {
                     updated_matching_order_bn.push_back(edge.first);
                 } else {
                     updated_matching_order_bn.push_back(edge.second);
                 }
             }
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else if (graph_size[selected_graph] == 2) {
             if (query_graph->getVertexDegree(connected_components[selected_graph][0]) < query_graph->getVertexDegree(connected_components[selected_graph][1])) {
                 std::swap(connected_components[selected_graph][0], connected_components[selected_graph][1]);
             }

             updated_matching_order.insert(updated_matching_order.end(), connected_components[selected_graph].begin(), connected_components[selected_graph].end());
             // Insert the first vertex.
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());

             // Insert the second vertex.
             updated_matching_order_bn.push_back(connected_components[selected_graph][0]);
             updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
             updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
         } else {
             auto graph = graphs[selected_graph];
             std::vector<bool> vertex_visited(matching_orders[selected_graph].size(), false);
             for (auto u: matching_orders[selected_graph]) {
                 uint32_t u_nbrs_count;
                 auto u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);

                 for (uint32_t j = 0; j < u_nbrs_count; ++j) {
                     uint32_t uu = u_nbrs[j];
                     if (vertex_visited[uu]) {
                         updated_matching_order_bn.push_back(connected_components[selected_graph][uu]);
                         updated_matching_order_edge_type.push_back(RelationEdgeType::REGULAR);
                     }
                 }

                 updated_matching_order.push_back(connected_components[selected_graph][u]);
                 updated_matching_order_bn_offset.push_back(updated_matching_order_bn.size());
                 vertex_visited[u] = true;
             }
         }
     }

    updated_matching_order_view_mapping.resize(updated_matching_order_bn.size());

     // Release graph.
     for (auto graph : graphs) {
         delete graph;
     }
}

void OrderManager::generate_matching_order_with_RI(const Graph *graph, std::vector<uint32_t> &matching_order) {
    uint32_t n = graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    // Select the vertex with the maximum degree as the start vertex.
    uint32_t selected_vertex = 0;
    uint32_t selected_vertex_selectivity = graph->getVertexDegree(selected_vertex);

    for (ui u = 1; u < n; ++u) {
        uint32_t u_selectivity = graph->getVertexDegree(u);
        if (u_selectivity > selected_vertex_selectivity) {
            selected_vertex = u;
            selected_vertex_selectivity = u_selectivity;
        }
    }

    matching_order.push_back(selected_vertex);
    visited[selected_vertex] = true;

    // Order vertices.
    std::vector<uint32_t> tie_vertices;
    std::vector<uint32_t> temp;

    for (uint32_t i = 1; i < n; ++i) {
        // Select the vertices with the maximum number of backward neighbors.
        selected_vertex_selectivity = 0;
        for (uint32_t u = 0; u < n; ++u) {
            if (!visited[u]) {
                // Compute the number of backward neighbors of u.
                uint32_t u_selectivity = 0;
                for (auto uu : matching_order) {
                    if (graph->checkEdgeExistence(u, uu)) {
                        u_selectivity += 1;
                    }
                }

                // Update the vertices under consideration.
                if (u_selectivity > selected_vertex_selectivity) {
                    selected_vertex_selectivity = u_selectivity;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                } else if (u_selectivity == selected_vertex_selectivity) {
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t count = 0;
            std::vector<uint32_t> u_fn;
            for (auto u : temp) {
                // Compute the number of vertices in the matching order that has at least one vertex
                // not in the matching order && connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count);
                for (uint32_t j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                uint32_t cur_count = 0;
                for (auto uu : matching_order) {
                    uint32_t uun_count;
                    auto uun = graph->getVertexNeighbors(uu, uun_count);
                    uint32_t common_neighbor_count = 0;
                    ComputeSetIntersection::ComputeCandidates(uun, uun_count, u_fn.data(), (uint32_t)u_fn.size(), common_neighbor_count);
                    if (common_neighbor_count > 0) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t count = 0;
            std::vector<uint32_t> u_fn;
            for (auto u : temp) {
                // Compute the number of vertices not in the matching order && not the neighbor of vertices in the
                // matching order, but is connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t un_count;
                auto un = graph->getVertexNeighbors(u, un_count);
                for (uint32_t j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                uint32_t cur_count = 0;
                for (auto uu : u_fn) {
                    bool valid = true;

                    for (auto uuu : matching_order) {
                        if (graph->checkEdgeExistence(uu, uuu)) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        matching_order.push_back(tie_vertices[0]);
        visited[tie_vertices[0]] = true;
        tie_vertices.clear();
        temp.clear();
    }
}

void OrderManager::create_orders(const Graph *query_graph) {
    // 1. Detect automorphisms.
    detect_automorphism_edges(query_graph);

    // 2. Create orders and map each edge to them.
    uint32_t count = 0;
    for (auto& automorphism : automorphism_edges_) {
        Edge edge = automorphism.front();
        orders_.emplace_back(OrdersPerEdge());

        create_indexing_order_for_reduced_graph(query_graph, edge, orders_.back());
        create_matching_order_for_reduced_graph(query_graph, edge, orders_.back());

        LabelTriple label_triple = { query_graph->getVertexLabel(edge.first),
                                     query_graph->getEdgeLabelByVertex(edge.first, edge.second), query_graph->getVertexLabel(edge.second)};
        {
            auto it = label_edge_mapping_.find(label_triple);
            if (it == label_edge_mapping_.end()) {
                auto temp_it = label_edge_mapping_.emplace(label_triple, std::vector<Edge>());
                it = temp_it.first;
            }

            for (auto e: automorphism) {
                edge_orders_mapping_.insert({e, orders_.size() - 1});
                it->second.push_back(e);
            }
        }
        {
            auto it = label_automorphism_mapping_.find(label_triple);
            if (it == label_automorphism_mapping_.end()) {
                auto temp_it = label_automorphism_mapping_.emplace(label_triple, std::vector<uint32_t>());
                it = temp_it.first;
            }
            it->second.push_back(orders_.size() - 1);
            count += 1;
        }
    }

//    printf("\n%zu\n", label_automorphism_mapping_.size());
//    printf("\n%zu\n", label_edge_mapping_.size());
//    printf("\n%u\n", count);
//    count = 0;
//    for (auto& x : label_edge_mapping_) {
//        count += x.second.size();
//    }
//    printf("\n%u\n", count);
//    count = 0;
//    for (auto& x : label_automorphism_mapping_) {
//        count += x.second.size();
//    }
//    printf("\n%u\n", count);
//
//    for (auto& kv : label_edge_mapping_) {
//        auto &test = kv.second;
//        printf("%u, %u, %u: ", kv.first.src_label_, kv.first.edge_label_, kv.first.dst_label_);
//        for (auto &x: test) {
//            printf("(%u, %u)\n", x.first, x.second);
//        }
//        printf("\n");
//    }
//
//    printf("%u, %u\n", query_graph->getEdgeLabelByVertex(4, 5), query_graph->getEdgeLabelByVertex(5, 4));
//    exit(-1);
}

void OrderManager::print_info() {
    printf("Order Info:\n");
    printf("The number of automorphisms: %zu\n", automorphisms_.size());

    for (auto& automorphism : automorphisms_) {
        for (auto u : automorphism) {
            printf("%u ", u);
        }
        printf("\n");
    }

    printf("-----\n");
    printf("The number of disjoint set: %zu\n", automorphism_edges_.size());
    for (uint32_t i = 0; i < automorphism_edges_.size(); ++i) {
        printf("Edge set: ");
        for (auto e : automorphism_edges_[i]) {
            printf("(%u, %u), ", e.first, e.second);
        }
        printf("\n");
        printf("Index order: ");
        for (auto u : orders_[i].indexing_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn offset: ");
        for (auto u : orders_[i].indexing_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Index order bn: ");
        for (auto u : orders_[i].indexing_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order: ");
        for (auto u : orders_[i].matching_order_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn offset: ");
        for (auto u : orders_[i].matching_order_bn_offset_) {
            printf("%u ", u);
        }
        printf("\n");

        printf("Matching order bn: ");
        for (auto u : orders_[i].matching_order_bn_) {
            printf("%u ", u);
        }
        printf("\n");
        printf("-----\n");
    }
}
