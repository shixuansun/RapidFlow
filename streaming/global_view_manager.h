//
// Created by sunsx on 28/05/21.
//

#ifndef RAPIDMATCH_GLOBAL_VIEW_MANAGER_H
#define RAPIDMATCH_GLOBAL_VIEW_MANAGER_H

#include "graph/graph.h"
#include "utility/relation/global_edge_view.h"
#include "utility/relation/encoded_edge_view.h"
#include "streaming_type.h"

class GlobalViewManager {
private:
    const Graph* query_graph_;
    std::vector<GlobalEdgeView> views_;
    spp::sparse_hash_map<LabelTriple, MappedViews> label_view_mapping_;
    spp::sparse_hash_map<Edge, uint32_t> edge_view_mapping_;
    uint32_t valid_view_count_;

    // NLF
    spp::sparse_hash_map<std::pair<LabelID, LabelID>, uint32_t> nlf_mapping_;
    std::vector<std::vector<uint32_t>> data_nlf_;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> query_nlf_;

    // Maintain the candidates for non-leaf query vertices.
    std::vector<spp::sparse_hash_map<uint32_t, uint32_t>> candidates_;
    std::vector<EncodedEdgeView> nlf_views_;
    spp::sparse_hash_map<Edge, uint32_t> edge_nlf_view_mapping_;
    std::vector<std::pair<uint32_t, uint32_t>> updated_candidates_;
    uint32_t valid_nlf_view_count_;

    // Encode the candidates.
    std::vector<std::vector<uint32_t>> encoded_candidates_;

public:
    GlobalViewManager() {}
    ~GlobalViewManager() { release(); }

    void initialize(const Graph *query_graph);
    void release();

    void create_views(const Graph *query_graph, const Graph *data_graph);

    MappedViews* get_mapped_views(LabelTriple label_triple);

    GlobalEdgeView* get_view(uint32_t view_id);

    GlobalEdgeView* get_view(Edge edge);

    EncodedEdgeView* get_nlf_view(Edge edge);

    uint32_t get_view_id(Edge edge);

    uint32_t get_encoded_candidate_set_size(uint32_t u);

    void update_view(char op, uint32_t view_id, Edge e, LabelTriple label_triple);

    void update_nlf_view(char op, Edge e, LabelTriple label_triple);

    void calculate_memory_usage();

    void print_view_info();

    inline bool nlf_check(uint32_t u, uint32_t v) {
        auto &u_nlf = query_nlf_[u];
        auto &v_nlf = data_nlf_[v];
        for (auto &kv : u_nlf) {
            if (v_nlf[kv.first] < kv.second)
                return false;
        }
        return true;
    }

    inline bool candidate_check(uint32_t u, uint32_t v) {
        return candidates_[u].contains(v);
    }

    spp::sparse_hash_map<uint32_t, uint32_t> *get_candidate_set(uint32_t u) {
        return &(candidates_[u]);
    }

    inline uint32_t get_id(uint32_t u, uint32_t encoded_id) {
        return encoded_candidates_[u][encoded_id];
    }

    inline uint32_t get_encoded_id(uint32_t u, uint32_t v) {
        return candidates_[u][v];
    }
};


#endif //RAPIDMATCH_GLOBAL_VIEW_MANAGER_H
