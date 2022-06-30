//
// Created by sunsx on 29/05/21.
//

#ifndef RAPIDMATCH_LOCAL_VIEW_MANAGER_H
#define RAPIDMATCH_LOCAL_VIEW_MANAGER_H

#include "streaming_type.h"
#include "global_view_manager.h"
#include "relation/local_edge_view.h"

class LocalViewManager {
public:
    uint64_t build_visited_neighbor_count_;
    uint64_t generate_visited_neighbor_count_;
    uint64_t first_vertex_neighbor_;
private:
    // Store the vertex id.
    std::vector<std::vector<uint32_t>> candidates_store_;
    std::vector<LocalEdgeView> views_;
    std::vector<uint32_t> buffer_pool_;
    spp::sparse_hash_map<Edge, uint32_t> edge_view_mapping_;

    std::vector<uint32_t> flag_array_;
    std::vector<uint32_t> updated_flag_;
    std::vector<uint32_t> si_buffer_;
    uint32_t updated_count_;

    std::vector<uint32_t*> candidate_set_pointer_;
    std::vector<uint32_t> candidate_set_size_;

private:
    bool optimized_generate_local_candidates(const Graph *query_graph, OrdersPerEdge &orders,
                                             GlobalViewManager &gvm, Edge exclude_data_edge);

    bool optimized_generate_local_candidates_v2(const Graph *query_graph, OrdersPerEdge &orders,
                                                GlobalViewManager &gvm, Edge exclude_data_edge);

    void optimized_build_local_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm);

    bool prune_local_candidates(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);

    void set_adjacent_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                             uint32_t encoded_v0, uint32_t encoded_v1);

    void set_update_candidates_flag(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn,
                                    uint32_t encoded_v0, uint32_t encoded_v1);

    uint32_t select_bn_with_minimum_degree_sum(GlobalViewManager &gvm, uint32_t u, std::vector<uint32_t> &bn);


public:
    LocalViewManager() {}
    ~LocalViewManager() { release(); }

    void initialize(const Graph *query_graph, const Graph *data_graph);
    void release();

    bool create_view(const Graph *query_graph, OrdersPerEdge &orders, GlobalViewManager &gvm, Edge data_edge);

    void destroy_view();

    LocalEdgeView* get_view(uint32_t view_id);
    LocalEdgeView* get_view(Edge query_edge);
    uint32_t get_view_id(Edge query_edge);
    uint32_t* get_candidate_set(uint32_t u);
    uint32_t get_candidate_set_size(uint32_t u);
};


#endif //RAPIDMATCH_LOCAL_VIEW_MANAGER_H
