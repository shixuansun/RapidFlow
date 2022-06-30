//
// Created by sunsx on 31/05/21.
//

#include <chrono>
#include "search_engine.h"
#include "computesetintersection.h"
#include "streaming_config.h"

void SearchEngine::initialize(const Graph *query_graph, const Graph *data_graph) {
    uint32_t n = query_graph->getVerticesCount();
    uint32_t N = data_graph->getVerticesCount();

    embedding_.resize(n);
    encoded_embedding_.resize(n);
    local_idx_.resize(n);
    local_candidates_store_.resize(n);
    encoded_local_candidates_store_.resize(n);

    local_candidates_buffer1_.resize(n);
    for (auto& buffer : local_candidates_buffer1_) {
        buffer.resize(N);
    }

    local_candidates_buffer2_.resize(n);
    for (auto& buffer : local_candidates_buffer2_) {
        buffer.resize(N);
    }

    visited_ = new bool [data_graph->getVerticesCount()];
    std::fill(visited_, visited_ + data_graph->getVerticesCount(), false);

    reset_performance_counters();
}

void SearchEngine::release() {
    embedding_.clear();
    encoded_embedding_.clear();
    local_idx_.clear();
    local_candidates_store_.clear();
    encoded_local_candidates_store_.clear();
    local_candidates_buffer1_.clear();
    local_candidates_buffer2_.clear();
    delete []visited_;
}



uint64_t SearchEngine::search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                               GlobalViewManager &gvm) {
    uint64_t result_count = 0;

    auto& order = orders.matching_order_;
    auto& bn_offset = orders.matching_order_bn_offset_;
    auto& bn = orders.matching_order_bn_;
    auto& view_mapping = orders.matching_order_view_mappings_;

    uint32_t target_depth = order.size();
    uint32_t start_vertex = order[0];
    uint32_t last_vertex = order[target_depth - 1];

    Edge data_edge = {lvm.get_candidate_set(orders.indexing_order_[0])[0], lvm.get_candidate_set(orders.indexing_order_[1])[0]};
    embedding_[orders.indexing_order_[0]] = data_edge.first;
    embedding_[orders.indexing_order_[1]] = data_edge.second;

    visited_[data_edge.first] = true;
    visited_[data_edge.second] = true;

    uint32_t* seeds;
    local_idx_[0].first = 0;
    if (query_graph->getVertexDegree(start_vertex) == 1) {
        local_idx_[0].second = compute_local_candidates_for_reduced_query(query_graph, 0, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);
        seeds = local_candidates_store_[0];
    }
    else {
        seeds = lvm.get_candidate_set(start_vertex);
        local_idx_[0].second = lvm.get_candidate_set_size(start_vertex);
    }



    if (target_depth == 1) {
#if EXECUTION_MODE == 1
        for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; ++local_idx_[0].first) {
            uint32_t v = (*seeds)[local_idx_[0].first];
            embedding_[start_vertex] = v;
        }
#endif
        visited_[data_edge.first] = false;
        visited_[data_edge.second] = false;
        return local_idx_[0].second;
    }


#ifdef ENABLE_PERFORMANCE_COUNTERS
    std::vector<uint64_t> result_count_vec(target_depth, 0);
#endif

    for (local_idx_[0].first = 0; local_idx_[0].first < local_idx_[0].second; local_idx_[0].first++) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
        partial_result_count_ += 1;
        result_count_vec[0] = result_count;
#endif

        uint32_t seed = seeds[local_idx_[0].first];
        if (visited_[seed])
            continue;

        embedding_[start_vertex] = seed;
        encoded_embedding_[start_vertex] = local_idx_[0].first;

        visited_[seed] = true;

        uint32_t current_depth = 1;
        local_idx_[current_depth].first = 0;
        local_idx_[current_depth].second = compute_local_candidates_for_reduced_query(query_graph, current_depth, order,
                                                                                      bn_offset, bn,
                                                                                      view_mapping, lvm, gvm);
        if (target_depth == 2) {
#if EXECUTION_MODE == 1
            for (;local_idx_[current_depth].first < local_idx_[current_depth].second; local_idx_[current_depth].first++) {
                uint32_t vv = local_candidates_store_[current_depth][local_idx_[current_depth].first];
                embedding_[last_vertex] = vv;
            }
#endif
            result_count += local_idx_[current_depth].second;
            visited_[seed] = false;

            if (result_count >= target_number)
                return result_count;
        } else {
            while (true) {
                while (local_idx_[current_depth].first < local_idx_[current_depth].second) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                    partial_result_count_ += 1;
                    result_count_vec[current_depth] = result_count;
#endif
                    uint32_t u = order[current_depth];
                    uint32_t encoded_v = encoded_local_candidates_store_[current_depth][local_idx_[current_depth].first];
                    uint32_t v = local_candidates_store_[current_depth][local_idx_[current_depth].first++];
                    encoded_embedding_[u] = encoded_v;
                    embedding_[u] = v;
                    visited_[v] = true;

                    uint32_t next_depth = current_depth + 1;
                    local_idx_[next_depth].first = 0;
                    local_idx_[next_depth].second = compute_local_candidates_for_reduced_query(query_graph, next_depth,
                                                                                               order,
                                                                                               bn_offset, bn,
                                                                                               view_mapping, lvm,
                                                                                               gvm);
                    if (local_idx_[next_depth].second == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
                        lc_empty_count_ += 1;
                        invalid_partial_result_count_ += 1;
#endif
                        visited_[v] = false;
                    } else if (next_depth == target_depth - 1) {
                        result_count += local_idx_[next_depth].second;
                        if (result_count >= target_number) {
                            for (int x = 0; x < target_depth; ++x) {
                                visited_[embedding_[order[x]]] = false;
                            }
                            visited_[data_edge.first] = false;
                            visited_[data_edge.second] = false;
                            return result_count;
                        }
////test

//                        if ((data_edge.first == 303387 && data_edge.second == 303399) || (data_edge.first == 303399 && data_edge.second == 303387)) {
//                            for (; local_idx_[next_depth].first <
//                                   local_idx_[next_depth].second; local_idx_[next_depth].first++) {
//                                uint32_t vv = local_candidates_store_[next_depth][local_idx_[next_depth].first];
//                                embedding_[last_vertex] = vv;
//                                result_count += 1;
//
//                                embedding_[orders.indexing_order_[0]] = lvm.get_candidate_set(
//                                        orders.indexing_order_[0])[0];
//                                embedding_[orders.indexing_order_[1]] = lvm.get_candidate_set(
//                                        orders.indexing_order_[1])[0];
//
//                                for (int x = 0; x < target_depth + 2; ++x) {
//                                    printf("(%u, %u), ", orders.indexing_order_[x],
//                                           embedding_[orders.indexing_order_[x]]);
//                                }
//                                printf("\n");
//                                if (result_count >= target_number) {
//                                    for (int x = 0; x < target_depth; ++x) {
//                                        visited_[embedding_[order[x]]] = false;
//                                    }
//                                    visited_[data_edge.first] = false;
//                                    visited_[data_edge.second] = false;
//                                    return result_count;
//                                }
//                            }
//                        }
//                            exit(-1);

////test
                        visited_[v] = false;

                        if(g_exit) {
                            return result_count;
                        }
                    } else {
                        current_depth += 1;
                    }
                }

                if (g_exit) {
                    return result_count;
                }

                current_depth -= 1;
                visited_[embedding_[order[current_depth]]] = false;
#ifdef ENABLE_PERFORMANCE_COUNTERS
                if (result_count == result_count_vec[current_depth]) {
                    invalid_partial_result_count_ += 1;
                }
#endif
                if (current_depth < 1) {
                    break;
                }
            }
        }
    }

    visited_[data_edge.first] = false;
    visited_[data_edge.second] = false;
    return result_count;
}

uint32_t SearchEngine::compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                                  std::vector<uint32_t> &order,
                                                                  std::vector<uint32_t> &bn_offset,
                                                                  std::vector<uint32_t> &bn,
                                                                  std::vector<uint32_t> &view_mapping,
                                                                  LocalViewManager &lvm, GlobalViewManager &gvm) {
    uint32_t bn_begin = bn_offset[depth];
    uint32_t bn_end = bn_offset[depth + 1];
    bool one_bn = bn_end - bn_begin == 1;
    uint32_t *lc1 = nullptr;
    uint32_t lc_count1 = 0;
    uint32_t *lc2 = nullptr;
    uint32_t lc_count2 = 0;

    if (bn_end - bn_begin == 0) {
        // Has no backward neighbors.
        lc1 = lvm.get_candidate_set(order[depth]);
        lc_count1 = lvm.get_candidate_set_size(order[depth]);
        lc2 = local_candidates_buffer1_[depth].data();

        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                local_candidates_buffer2_[depth][lc_count2] = i;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        uint32_t u = bn[bn_begin];
        uint32_t v = embedding_[u];
        uint32_t encoded_v = encoded_embedding_[u];

        uint32_t view_id = view_mapping[bn_begin];

        // If the vertex degree is 1 (i.e., a leaf), retrieve the neighbor from global view.
        lc1 = query_graph->getVertexDegree(order[depth]) == 1 ? gvm.get_view(view_id)->get_neighbor(v, lc_count1)
                : lvm.get_view(view_id)->get_neighbors(encoded_v, lc_count1);
//        lc1 = lvm.get_view(view_id)->get_neighbors(v, lc_count1);

        lc2 = local_candidates_buffer1_[depth].data();

        if (lc_count1 == 0) {
            goto EXIT;
        }

        // More than one backward neighbor.
        if (bn_begin + 1 < bn_end) {
            bn_begin += 1;
            u = bn[bn_begin];
            encoded_v = encoded_embedding_[u];
            view_id = view_mapping[bn_begin];

            lc2 = lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);

            uint32_t temp_count;
            uint32_t *temp_buffer = local_candidates_buffer1_[depth].data();

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

            if (temp_count == 0) {
                lc_count1 = 0;
                goto EXIT;
            }

            lc1 = temp_buffer;
            lc_count1 = temp_count;
            temp_buffer = local_candidates_buffer2_[depth].data();

            for (bn_begin += 1; bn_begin < bn_end; ++bn_begin) {
                u = bn[bn_begin];
                encoded_v = encoded_embedding_[u];
                view_id = view_mapping[bn_begin];

                lc2 =  lvm.get_view(view_id)->get_neighbors(encoded_v,lc_count2);
                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, temp_buffer, temp_count);

                if (temp_count == 0) return 0;
                std::swap(temp_buffer, lc1);
                std::swap(temp_count, lc_count1);
            }

            lc2 = temp_buffer;
        }
    }

EXIT:
    if (lc_count1 == 0) {
#ifdef ENABLE_PERFORMANCE_COUNTERS
       si_empty_count_ += 1;
#endif
        return 0;
    }

    if (query_graph->getVertexDegree(order[depth]) == 1) {
        lc_count2 = 0;

        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t v = lc1[i];

            if (!visited_[v]) {
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }
        encoded_local_candidates_store_[depth] = local_candidates_buffer2_[depth].data();
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
    else {
        lc_count2 = 0;
        uint32_t* candidate_set = lvm.get_candidate_set(order[depth]);
        uint32_t* encoded_buffer = one_bn ? local_candidates_buffer2_[depth].data() : lc1;
        for (uint32_t i = 0; i < lc_count1; ++i) {
            uint32_t encoded_v = lc1[i];
            uint32_t v = candidate_set[encoded_v];

            if (!visited_[v]) {
                encoded_buffer[lc_count2] = encoded_v;
                lc2[lc_count2++] = v;
            }
#ifdef ENABLE_PERFORMANCE_COUNTERS
            else {
                iso_conflict_count_ += 1;
            }
#endif
        }

        encoded_local_candidates_store_[depth] = encoded_buffer;
        local_candidates_store_[depth] = lc2;
        return lc_count2;
    }
}


