//
// Created by sunsx on 02/06/21.
//

#ifndef RAPIDMATCH_STREAMING_ENGINE_H
#define RAPIDMATCH_STREAMING_ENGINE_H

#include "order_manager.h"
#include "search_engine.h"

class StreamingEngine {
public:
    // performance counters.
    bool is_relevant_;
    bool is_searched_;

    uint64_t query_time_ = 0;
    uint64_t global_view_initialize_time_ = 0;
    uint64_t order_generation_time_ = 0;

    // The average time of executing 5 times.
    uint64_t global_view_update_time_ = 0;

    // The average time of executing 5 times.
    uint64_t local_view_update_time_ = 0;

    // The number of edges processed.
    uint64_t edge_process_count_ = 0;

    // The number of updates relevant to the query.
    uint64_t relevant_update_count_ = 0;

    // The number of updates triggering the search.
    uint64_t search_count_ = 0;

    // The number of updates leading to new results.
    uint64_t positive_count_ = 0;

    uint64_t result_count_ = 0;
    uint64_t invalid_partial_result_count_ = 0;
    uint64_t partial_result_count_ = 0;
    uint64_t iso_conflict_count_ = 0;
    uint64_t si_empty_count_ = 0;
    uint64_t lc_empty_count_ = 0;

    uint64_t non_search_generate_neighbor_count_;
    uint64_t search_generate_neighbor_count_;
    uint64_t non_search_build_neighbor_count_;
    uint64_t search_build_neighbor_count_;
    uint64_t direct_rejection_count_;
    uint64_t first_indexing_vertex_;
    void reset_performance_counters() {
        edge_process_count_ = 0;
        relevant_update_count_ = 0;
        search_count_ = 0;
        positive_count_ = 0;

        result_count_ = 0;
        invalid_partial_result_count_ = 0;
        partial_result_count_ = 0;
        iso_conflict_count_ = 0;
        si_empty_count_ = 0;
        lc_empty_count_ = 0;
    }
private:
    OrderManager om_;
    GlobalViewManager gvm_;
    LocalViewManager lvm_;
    SearchEngine sm_;

public:
    void initialize(const Graph *query_graph, const Graph *data_graph, uint64_t target_embedding_num);
    void preprocess(const Graph* query_graph, const Graph* data_graph);
    uint64_t execute(const Graph *query_graph, const Update &update, bool enable_local_view, bool enable_search);
    void release();

    void evaluate_view_update(const Graph *query_graph, const Graph *data_graph,
                              const std::vector<Update> &stream);
    void evaluate_search(const Graph *query_graph, const Graph *data_graph, const std::vector<Update> &stream);

    void print_metrics();
};

#endif //RAPIDMATCH_STREAMING_ENGINE_H
