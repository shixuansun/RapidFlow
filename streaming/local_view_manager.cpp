//
// Created by sunsx on 29/05/21.
//

#include <chrono>
#include <cmath>
#include "local_view_manager.h"
#include "streaming_config.h"
#include "computesetintersection.h"
void LocalViewManager::initialize(const Graph *query_graph, const Graph *data_graph) {
    candidates_store_.resize(query_graph->getVerticesCount());
    for (auto& candidates_set : candidates_store_) {
        candidates_set.reserve(1024);
    }
    views_.resize(query_graph->getEdgesCount() * 2);
    buffer_pool_.reserve(1024 * 1024);
    edge_view_mapping_.reserve(256);
    flag_array_.resize(data_graph->getVerticesCount(), 0);
    updated_flag_.resize(data_graph->getVerticesCount());
    si_buffer_.resize(data_graph->getVerticesCount());
    candidate_set_size_.resize(query_graph->getVerticesCount());
    candidate_set_pointer_.resize(query_graph->getVerticesCount());
}

void LocalViewManager::release() {
    candidates_store_.clear();
    views_.clear();
    buffer_pool_.clear();
    edge_view_mapping_.clear();
    flag_array_.clear();
    updated_flag_.clear();
    candidate_set_size_.clear();
    candidate_set_pointer_.clear();
    si_buffer_.clear();
}

void LocalViewManager::destroy_view() {
    for (auto& view : views_) {
        view.clear();
    }

    for (auto& candidate_set : candidates_store_) {
        candidate_set.clear();
    }

    buffer_pool_.clear();
    edge_view_mapping_.clear();
}

LocalEdgeView *LocalViewManager::get_view(uint32_t view_id) {
    assert(view_id < views_.size());
    return &views_[view_id];
}

LocalEdgeView *LocalViewManager::get_view(Edge query_edge) {
    auto it = edge_view_mapping_.find(query_edge);
    if (it != edge_view_mapping_.end()) {
        return &views_[it->second];
    }
    return nullptr;
}

uint32_t LocalViewManager::get_view_id(Edge query_edge) {
    auto it = edge_view_mapping_.find(query_edge);
    if (it != edge_view_mapping_.end()) {
        return it->second;
    }
    // Local view does not exist.
    return edge_view_mapping_.size();
}

uint32_t *LocalViewManager::get_candidate_set(uint32_t u) {
    return candidates_store_[u].data();
}

uint32_t LocalViewManager::get_candidate_set_size(uint32_t u) {
    return candidates_store_[u].size();
}

bool
LocalViewManager::create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge) {
//    auto start = std::chrono::high_resolution_clock::now();
    build_visited_neighbor_count_  = 0;
    generate_visited_neighbor_count_ = 0;
    first_vertex_neighbor_ = 0;
//    if (orders.indexing_order_[0] == 5 && orders.indexing_order_[1] == 1
//        && data_edge.first == 2765688 && data_edge.second == 2183133) {
//        printf("Hello\n");
//    }

    bool is_valid = optimized_generate_local_candidates_v2(query_graph, orders, gvm, data_edge);

    if (!is_valid)
        return false;
//    auto end = std::chrono::high_resolution_clock::now();
//    auto generate_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
//    printf("Generate time (seconds): %f\n", NANOSECTOSEC(generate_time));
//    printf("(%d, %d) (%d, %d)\n", orders.indexing_order_[0], orders.indexing_order_[1], data_edge.first, data_edge.second);
//     start = std::chrono::high_resolution_clock::now();
     optimized_build_local_view(query_graph, orders, gvm);

    // build_local_view(query_graph, orders, gvm);
//     end = std::chrono::high_resolution_clock::now();
//     auto build_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
//     printf("Build time (seconds): %f\n", NANOSECTOSEC(build_time));
//     printf("Build visited neighbor count: %zu\n", build_visited_neighbor_count_);
//    uint32_t count = 0;
//    for (auto u : orders.indexing_order_) {
//        count += candidates_store_[u].size();
//        printf("(%u, %zu) ", u, candidates_store_[u].size());
//    }
//    printf("\nTotal %u\n", count);

    return true;
}

bool LocalViewManager::optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders,
                                                           GlobalViewManager &gvm, Edge exclude_data_edge) {
    // Generate the candidate set for each vertex along the indexing order.
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;

    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;

    if (!gvm.nlf_check(u0, v0)
        || !gvm.nlf_check(u1, v1)) {
        return false;
    }

    if (query_graph->getVertexDegree(u0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(u0, v0);
        candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        candidates_store_[u0].push_back(v0);
    }

    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(u1, v1);
        candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        candidates_store_[u1].push_back(v1);
    }

    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u = indexing_order[i];

        // Skip leaf node.
        if (query_graph->getVertexDegree(u) == 1)
            continue;

        uint32_t begin = indexing_order_bn_offset[i];
        uint32_t end = indexing_order_bn_offset[i + 1];
        updated_count_ = 0;

        uint32_t flag_value = 0;
        for (uint32_t j = begin; j < end; ++j) {
            uint32_t uu = indexing_order_bn[j];
            auto gv = gvm.get_nlf_view({uu, u});
            auto reverse_gv = gvm.get_nlf_view({u, uu});

            bool flag = false;
            if (end - begin >= 2) {
                flag = (indexing_order_bn[begin] == u0 && indexing_order_bn[begin + 1] == u1);
            }
            for (auto vv : candidates_store_[uu]) {
                uint32_t vv_nbrs_count;
                auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                generate_visited_neighbor_count_ += vv_nbrs_count;
                if ((uu == indexing_order[0] || uu == indexing_order[1]) && flag)
                    first_vertex_neighbor_ += vv_nbrs_count;

                // If it is the first bn or the cost of binary search is lower than that of scan, then use the first approach.
                if ((j == begin) || vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                    for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                        uint32_t v = vv_nbrs[k];

                        if (flag_array_[v] == flag_value) {
                            flag_array_[v] += 1;
                            if (flag_value == 0) {
                                updated_flag_[updated_count_++] = v;
                            }
                        }
                    }
                }
                else {
                    for (uint32_t k = 0; k < updated_count_; ++k) {
                        uint32_t v = updated_flag_[k];
                        if (flag_array_[v] == flag_value) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    flag_array_[v] += 1;
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    flag_array_[v] += 1;
                                }
                            }
                        }
                    }
                }
            }

            flag_value += 1;
        }
        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(u, v0);
        }

        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(u, v1);
        }

        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = updated_flag_[j];
            if (flag_array_[v] == flag_value && v != local_encoded_v0 && v != local_encoded_v1) {
                candidates_store_[u].push_back(v);
            }

            flag_array_[v] = 0;
        }

        if (candidates_store_[u].empty())
            return false;
    }

    for (auto& candidate_set : candidates_store_) {
        std::sort(candidate_set.begin(), candidate_set.end());
    }

    return true;
}

void
LocalViewManager::optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm) {
    // Create views based on the candidate set along the matching order.
    auto& matching_order = orders.matching_order_;
    auto& matching_order_bn = orders.matching_order_bn_;
    auto& matching_order_bn_offset = orders.matching_order_bn_offset_;
    auto& matching_order_view_mappings = orders.matching_order_view_mappings_;
    auto& matching_order_relation_edge_type = orders.matching_order_edge_type_;

    uint32_t local_view_id = 0;
    for (uint32_t i = 0; i < matching_order.size(); ++i) {
        uint32_t u = matching_order[i];

        updated_count_ = 0;

        for (uint32_t j = matching_order_bn_offset[i]; j < matching_order_bn_offset[i + 1]; ++j) {
            uint32_t uu = matching_order_bn[j];
            RelationEdgeType type = matching_order_relation_edge_type[j];
            Edge query_edge = {uu, u};
            Edge reverse_query_edge = {u, uu};

            if (query_graph->getVertexDegree(u) == 1) {
                matching_order_view_mappings[j]= gvm.get_view_id(query_edge);
                continue;
            }

            matching_order_view_mappings[j] = local_view_id;
            edge_view_mapping_[query_edge] = local_view_id;

            auto &lv = views_[local_view_id++];

            if (type == RelationEdgeType::REGULAR) {
                // Set the flag array.
                if (updated_count_ == 0) {
                    for (uint32_t k = 0; k < candidates_store_[u].size(); ++k) {
                        flag_array_[candidates_store_[u][k]] = k + 1;
                        updated_flag_[updated_count_++] = candidates_store_[u][k];
                    }
                }

                // Generate the local view.
                auto gv = gvm.get_nlf_view(query_edge);
                auto reverse_gv = gvm.get_nlf_view(reverse_query_edge);

                for (auto vv : candidates_store_[uu]) {
                    uint32_t begin_pos = buffer_pool_.size();

                    uint32_t vv_nbrs_count;
                    auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);
                    build_visited_neighbor_count_ += vv_nbrs_count;

                    if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                        for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                            uint32_t v = vv_nbrs[k];
                            if (flag_array_[v] > 0) {
                                buffer_pool_.push_back(flag_array_[v] - 1);
                            }
                        }
                    }
                    else {
 //                       printf("B %u, %u, %d, %d\n", candidates_store_[u].size(), vv_nbrs_count, updated_count_);
                        for (auto v : candidates_store_[u]) {
                            uint32_t v_nbrs_count;
                            auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                            if (vv_nbrs_count < v_nbrs_count) {
                                auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                                if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                    buffer_pool_.push_back(flag_array_[v] - 1);
                                }
                            }
                            else {
                                auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                                if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                    buffer_pool_.push_back(flag_array_[v] - 1);
                                }
                            }
                        }
                    }

                    lv.trie_.emplace_back(begin_pos, buffer_pool_.size());
                    lv.cardinality_ += buffer_pool_.size() - begin_pos;
                }
            }
        }

        for (uint32_t j = 0; j < updated_count_; ++j) {
            uint32_t v = updated_flag_[j];
            flag_array_[v] = 0;
        }
    }

    // NOTE: Set the pointer in the local views at last, because the memory is dynamically allocated.
    for (uint32_t i = 0; i < local_view_id; ++i) {
        views_[i].data_ = buffer_pool_.data();
    }

    // Decode the local candidates.
    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            for (uint32_t i = 0; i < candidates_store_[u].size(); ++i) {
                candidates_store_[u][i] = gvm.get_id(u, candidates_store_[u][i]);
            }
        }
    }
}

bool LocalViewManager::optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                              GlobalViewManager &gvm, Edge exclude_data_edge) {
    // Generate the candidate set for each vertex along the indexing order.
    auto& indexing_order = orders.indexing_order_;
    auto& indexing_order_bn = orders.indexing_order_bn_;
    auto& indexing_order_bn_offset = orders.indexing_order_bn_offset_;

    uint32_t u0 = indexing_order[0];
    uint32_t u1 = indexing_order[1];
    uint32_t v0 = exclude_data_edge.first;
    uint32_t v1 = exclude_data_edge.second;

    if (!gvm.nlf_check(u0, v0)
        || !gvm.nlf_check(u1, v1)) {
        return false;
    }

    if (query_graph->getVertexDegree(u0) > 1) {
        uint32_t global_encoded_v0 = gvm.get_encoded_id(u0, v0);
        candidates_store_[u0].push_back(global_encoded_v0);
    }
    else {
        candidates_store_[u0].push_back(v0);
    }

    if (query_graph->getVertexDegree(u1) > 1) {
        uint32_t global_encoded_v1 = gvm.get_encoded_id(u1, v1);
        candidates_store_[u1].push_back(global_encoded_v1);
    }
    else {
        candidates_store_[u1].push_back(v1);
    }

    candidate_set_size_[u0] = 1;
    candidate_set_pointer_[u0] = candidates_store_[u0].data();
    candidate_set_size_[u1] = 1;
    candidate_set_pointer_[u1] = candidates_store_[u1].data();

    // First, compute the triangle.
    for (uint32_t i = 2; i < orders.triangle_end_index_; ++i) {
        // Skip the first two vertex.
        uint32_t u = indexing_order[i];

        auto gv0 = gvm.get_nlf_view({u0, u});
        uint32_t v0_nbr_count;
        auto v0_nbr = gv0->get_neighbor(candidate_set_pointer_[u0][0], v0_nbr_count);

        auto gv1 = gvm.get_nlf_view({u1, u});
        uint32_t v1_nbr_count;
        auto v1_nbr = gv1->get_neighbor(candidate_set_pointer_[u1][0], v1_nbr_count);

        uint32_t lc = 0;
        ComputeSetIntersection::ComputeCandidates(v0_nbr, v0_nbr_count, v1_nbr, v1_nbr_count, updated_flag_.data(), lc);

        if (lc == 0)
            return false;

        candidates_store_[u].insert(candidates_store_[u].end(), updated_flag_.begin(), updated_flag_.begin() + lc);
        candidate_set_pointer_[u] = candidates_store_[u].data();
        candidate_set_size_[u] = candidates_store_[u].size();
    }

    // Second, set the candidate for vertex adjacent to update.
    for (uint32_t i = orders.triangle_end_index_; i < orders.adjacent_end_index_; ++i) {
        uint32_t u = indexing_order[i];

        // Skip leaf node.
        if (query_graph->getVertexDegree(u) == 1)
            continue;

        uint32_t uu = indexing_order_bn[indexing_order_bn_offset[i]];
        auto gv = gvm.get_nlf_view({uu, u});
        candidate_set_pointer_[u] = gv->get_neighbor(candidate_set_pointer_[uu][0], candidate_set_size_[u]);

        if (candidate_set_size_[u] < 1024) {
            uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(u, v0)) {
                local_encoded_v0 = gvm.get_encoded_id(u, v0);
            }

            uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
            if (gvm.candidate_check(u, v1)) {
                local_encoded_v1 = gvm.get_encoded_id(u, v1);
            }

            for (uint32_t j = 0; j < candidate_set_size_[u]; ++j) {
                uint32_t v = candidate_set_pointer_[u][j];
                if (v != local_encoded_v0 && v != local_encoded_v1) {
                    candidates_store_[u].push_back(v);
                }
            }

            candidate_set_pointer_[u] = candidates_store_[u].data();
            candidate_set_size_[u] = candidates_store_[u].size();
        }

        if (candidate_set_size_[u] == 0)
            return false;
    }

        /*
         * Regenerate indexing order for adjacent vertex based on candidate set size.
         */
        std::vector<uint32_t> optimized_indexing_order(indexing_order.begin() + 2,
                                                       indexing_order.begin() + orders.adjacent_end_index_);
        std::sort(optimized_indexing_order.begin(), optimized_indexing_order.end(), [this](uint32_t a, uint32_t b) {
            return candidate_set_size_[a] < candidate_set_size_[b];
        });

        std::vector<bool> visited(query_graph->getVerticesCount(), false);
        std::vector<uint32_t> optimized_indexing_order_bn;
        std::vector<uint32_t> optimized_indexing_order_bn_offset;

        optimized_indexing_order_bn_offset.push_back(0);
        for (auto u : optimized_indexing_order) {
            uint32_t u_nbrs_count;
            auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
            for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                uint32_t uu = u_nbrs[i];
                if (visited[uu]) {
                    optimized_indexing_order_bn.push_back(uu);
                }
            }
            optimized_indexing_order_bn_offset.push_back(optimized_indexing_order_bn.size());
            visited[u] = true;
        }


    for (uint32_t i = 2; i < indexing_order.size(); ++i) {
        uint32_t u;
        std::vector<uint32_t> bn;
        if (i < orders.adjacent_end_index_) {
            uint32_t index = i - 2;
            u = optimized_indexing_order[index];
            uint32_t begin = optimized_indexing_order_bn_offset[index];
            uint32_t end = optimized_indexing_order_bn_offset[index + 1];
            bn.insert(bn.end(), optimized_indexing_order_bn.begin() + begin, optimized_indexing_order_bn.begin() + end);
        }
        else {
            u = indexing_order[i];
            uint32_t begin = indexing_order_bn_offset[i];
            uint32_t end = indexing_order_bn_offset[i + 1];
            bn.insert(bn.end(), indexing_order_bn.begin() + begin, indexing_order_bn.begin() + end);
        }
        if (query_graph->getVertexDegree(u) == 1 || bn.empty())
            continue;


        uint32_t local_encoded_v0 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v0)) {
            local_encoded_v0 = gvm.get_encoded_id(u, v0);
        }

        uint32_t local_encoded_v1 = std::numeric_limits<uint32_t>::max();
        if (gvm.candidate_check(u, v1)) {
            local_encoded_v1 = gvm.get_encoded_id(u, v1);
        }

        if (i < orders.adjacent_end_index_) {
            set_adjacent_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }
        else {
            set_update_candidates_flag(gvm, u, bn, local_encoded_v0, local_encoded_v1);
        }

        if (!prune_local_candidates(gvm, u, bn))
            return false;
    }

    for (auto u : indexing_order) {
        if (query_graph->getVertexDegree(u) > 1) {
            if (candidates_store_[u].empty()) {
                for (uint32_t i = 0; i < candidate_set_size_[u]; ++i) {
                    candidates_store_[u].push_back(candidate_set_pointer_[u][i]);
                }
            }
            std::sort(candidates_store_[u].begin(), candidates_store_[u].end());
            candidate_set_pointer_[u] = candidates_store_[u].data();
            candidate_set_size_[u] = candidates_store_[u].size();
        }
    }
    return true;
}

bool LocalViewManager::prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    uint32_t flag_value = 1;
    uint32_t current_valid_count = updated_count_;

    for (auto uu : bn) {
        auto gv = gvm.get_nlf_view({uu, u});
        auto reverse_gv = gvm.get_nlf_view({u, uu});

        uint32_t target_valid_count = current_valid_count;
        current_valid_count = 0;

        for (uint32_t j = 0; j < candidate_set_size_[uu]; ++j) {
            uint32_t vv = candidate_set_pointer_[uu][j];
            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            // If the cost of binary search is lower than that of scan, then use the first approach.
            if (vv_nbrs_count < 1024 || vv_nbrs_count < updated_count_ * 32) {
                for (uint32_t k = 0; k < vv_nbrs_count; ++k) {
                    uint32_t v = vv_nbrs[k];

                    if (flag_array_[v] == flag_value) {
                        flag_array_[v] += 1;
                        current_valid_count += 1;
                    }
                }
            }
            else {
                for (uint32_t k = 0; k < updated_count_; ++k) {
                    uint32_t v = updated_flag_[k];
                    if (flag_array_[v] == flag_value) {
                        uint32_t v_nbrs_count;
                        auto v_nbrs = reverse_gv->get_neighbor(v, v_nbrs_count);
                        if (vv_nbrs_count < v_nbrs_count) {
                            auto it = std::lower_bound(vv_nbrs, vv_nbrs + vv_nbrs_count, v);
                            if (it != vv_nbrs + vv_nbrs_count && *it == v) {
                                flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                        else {
                            auto it = std::lower_bound(v_nbrs, v_nbrs + v_nbrs_count, vv);
                            if (it != v_nbrs + v_nbrs_count && *it == vv) {
                                flag_array_[v] += 1;
                                current_valid_count += 1;
                            }
                        }
                    }
                }
            }

            if (current_valid_count == target_valid_count)
                break;
        }

        flag_value += 1;
    }

    bool push = candidates_store_[u].empty();
    uint32_t local_pos = 0;
    for (uint32_t j = 0; j < updated_count_; ++j) {
        uint32_t v = updated_flag_[j];
        if (flag_array_[v] == flag_value) {
            if (push)
                candidates_store_[u].push_back(v);
            else
                candidates_store_[u][local_pos++] = v;
        }

        flag_array_[v] = 0;
    }

    if (!push)
        candidates_store_[u].resize(local_pos);

    candidate_set_pointer_[u] = candidates_store_[u].data();
    candidate_set_size_[u] = candidates_store_[u].size();

    if (candidate_set_size_[u] == 0)
        return false;

    return true;
}

void
LocalViewManager::set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                      uint32_t encoded_v0, uint32_t encoded_v1) {
    bool is_flag_based = false;
    if (candidate_set_size_[u] < 1024) {
        is_flag_based = true;
    }
    else {
        for (auto uu : bn) {
            if (candidate_set_size_[u] < 32 * candidate_set_size_[uu]) {
                is_flag_based = true;
                break;
            }
        }
    }

    if (is_flag_based) {
        updated_count_ = 0;
        for (uint32_t i = 0; i < candidate_set_size_[u]; ++i) {
            uint32_t v = candidate_set_pointer_[u][i];
            if (v != encoded_v0 && v != encoded_v1) {
                flag_array_[v] = 1;
                updated_flag_[updated_count_++] = v;
            }
        }
    }
    else {
        uint32_t uu = select_bn_with_minimum_degree_sum(gvm, u, bn);
        auto gv = gvm.get_nlf_view({uu, u});

        updated_count_ = 0;
        for (uint32_t i = 0; i < candidate_set_size_[uu]; ++i) {
            uint32_t vv = candidate_set_pointer_[uu][i];

            uint32_t vv_nbrs_count;
            auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

            uint32_t lc = 0;
            ComputeSetIntersection::ComputeCandidates(candidate_set_pointer_[u], candidate_set_size_[u],
                                                      vv_nbrs, vv_nbrs_count, si_buffer_.data(), lc);


            for (uint32_t j = 0; j < lc; ++j) {
                uint32_t v = si_buffer_[j];

                if (flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                    flag_array_[v] = 1;
                    updated_flag_[updated_count_++] = v;
                }
            }
            if (candidate_set_size_[u] == updated_count_)
                break;
        }
    }
}

void LocalViewManager::set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                                  uint32_t encoded_v0, uint32_t encoded_v1) {

    uint32_t selected_bn = select_bn_with_minimum_degree_sum(gvm, u, bn);
    auto gv = gvm.get_nlf_view({selected_bn, u});
    updated_count_ = 0;
    for (uint32_t i = 0; i < candidate_set_size_[selected_bn]; ++i) {
        uint32_t vv = candidate_set_pointer_[selected_bn][i];
        uint32_t vv_nbrs_count;
        auto vv_nbrs = gv->get_neighbor(vv, vv_nbrs_count);

        for (uint32_t j = 0; j < vv_nbrs_count; ++j) {
            uint32_t v = vv_nbrs[j];
            if (flag_array_[v] == 0 && v != encoded_v0 && v != encoded_v1) {
                flag_array_[v] = 1;
                updated_flag_[updated_count_++] = v;
            }
        }
    }
}

uint32_t
LocalViewManager::select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn) {
    std::sort(bn.begin(), bn.end(), [this](const uint32_t a, const uint32_t b) -> bool {
        return candidate_set_size_[a] < candidate_set_size_[b];
    });

    uint64_t min_degree_sum = std::numeric_limits<uint64_t>::max();
    uint32_t selected_bn = 0;
    uint32_t selected_idx = 0;
    for (uint32_t i = 0; i < bn.size(); ++i) {
        uint32_t uu = bn[i];
        auto gv = gvm.get_nlf_view({uu, u});
        uint64_t local_degree_sum = 0;
        if (candidate_set_size_[uu] < min_degree_sum) {
            for (uint32_t j = 0; j < candidate_set_size_[uu]; ++j) {
                uint32_t vv = candidate_set_pointer_[uu][j];
                local_degree_sum += gv->get_neighbor_num(vv);
            }

            if (local_degree_sum < min_degree_sum) {
                min_degree_sum = local_degree_sum;
                selected_bn = uu;
                selected_idx = i;
            }
        }
    }

    bn.erase(bn.begin() + selected_idx);
    return selected_bn;
}