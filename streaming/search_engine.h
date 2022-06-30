//
// Created by sunsx on 31/05/21.
//

#ifndef RAPIDMATCH_SEARCH_ENGINE_H
#define RAPIDMATCH_SEARCH_ENGINE_H


#include "local_view_manager.h"

using spp::sparse_hash_set;

class SearchEngine {
public:
    uint64_t target_number = 1000;
public:
    // Performance counters
    uint64_t invalid_partial_result_count_;
    uint64_t partial_result_count_;
    uint64_t iso_conflict_count_;
    uint64_t si_empty_count_;
    uint64_t lc_empty_count_;

    void reset_performance_counters() {
        invalid_partial_result_count_ = 0;
        partial_result_count_ = 0;
        iso_conflict_count_ = 0;
        si_empty_count_ = 0;
        lc_empty_count_ = 0;
    }
private:
    std::vector<uint32_t*> local_candidates_store_;
    std::vector<uint32_t*> encoded_local_candidates_store_;

    std::vector<std::vector<uint32_t>> local_candidates_buffer1_;
    std::vector<std::vector<uint32_t>> local_candidates_buffer2_;
    std::vector<std::pair<uint32_t, uint32_t>> local_idx_;
    // spp::sparse_hash_set<uint32_t> visited_;
    bool* visited_;

    // Map u to v
    std::vector<uint32_t> embedding_;

    // Map u to local encoded id
    std::vector<uint32_t> encoded_embedding_;

private:
    uint32_t compute_local_candidates_for_reduced_query(const Graph *query_graph, uint32_t depth,
                                                        std::vector<uint32_t> &order,
                                                        std::vector<uint32_t> &bn_offset,
                                                        std::vector<uint32_t> &bn,
                                                        std::vector<uint32_t> &view_mapping,
                                                        LocalViewManager &lvm, GlobalViewManager &gvm);
public:
    SearchEngine() {}
    ~SearchEngine() {}

    void initialize(const Graph *query_graph, const Graph *data_graph);
    void release();

    uint64_t search_on_reduced_query(const Graph *query_graph, OrdersPerEdge &orders, LocalViewManager &lvm,
                                     GlobalViewManager &gvm);
};


#endif //RAPIDMATCH_SEARCH_ENGINE_H
