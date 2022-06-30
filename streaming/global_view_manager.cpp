//
// Created by sunsx on 28/05/21.
//

#include "global_view_manager.h"

void GlobalViewManager::initialize(const Graph *query_graph) {
    query_graph_ = query_graph;
    views_.resize(query_graph->getEdgesCount() * 2);
    valid_view_count_ = 0;

    uint32_t n = query_graph->getVerticesCount();
    query_nlf_.resize(n);
    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_nbrs_count;
        auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        spp::sparse_hash_map<std::pair<LabelID, LabelID>, uint32_t> frequency;

        for (uint32_t i = 0; i < u_nbrs_count; ++i) {
            uint32_t uu = u_nbrs[i];
            std::pair<LabelID, LabelID> key = {query_graph->getEdgeLabelByVertex(u, uu), query_graph->getVertexLabel(uu)};
            if (!nlf_mapping_.contains(key)) {
                nlf_mapping_[key] = nlf_mapping_.size();
                frequency[key] = 0;
            }

            frequency[key] += 1;
        }

        for (auto& kv : frequency) {
            query_nlf_[u].emplace_back(nlf_mapping_[kv.first], kv.second);
        }
    }

    candidates_.resize(n);
    nlf_views_.resize(query_graph->getEdgesCount() * 2);
    valid_nlf_view_count_ = 0;

    encoded_candidates_.resize(n);
}

void GlobalViewManager::release() {
    query_graph_ = nullptr;
    views_.clear();
    label_view_mapping_.clear();
    edge_view_mapping_.clear();
    valid_view_count_ = 0;
    data_nlf_.clear();
    query_nlf_.clear();
    nlf_mapping_.clear();
    candidates_.clear();
    nlf_views_.clear();
    edge_nlf_view_mapping_.clear();
    updated_candidates_.clear();
    encoded_candidates_.clear();
    valid_nlf_view_count_ = 0;
}

void GlobalViewManager::create_views(const Graph *query_graph, const Graph *data_graph) {
    // Initialize nlf.
    uint32_t N = data_graph->getVerticesCount();
    data_nlf_.resize(N);
    for (auto& nlf : data_nlf_) {
        nlf.resize(nlf_mapping_.size(), 0);
    }

    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        uint32_t u_nbrs_count;
        auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

        for (uint32_t i = 0; i < u_nbrs_count; ++i) {
            uint32_t uu = u_nbrs[i];
            auto u_label = query_graph->getVertexLabel(u);
            auto uu_label = query_graph->getVertexLabel(uu);

            uint32_t query_edge_label = query_graph->getEdgeLabelByLocalOffset(u, i);

            LabelTriple label_triple = {u_label, query_edge_label, uu_label};
            Edge e = {u, uu};

            // check whether the global view exists.
            auto view_it = label_view_mapping_.find(label_triple);
            if (view_it != label_view_mapping_.end()) {
                // If exist, then update the mapped view.
                view_it->second.first.push_back(e);
                edge_view_mapping_[e] = view_it->second.second;
            }
            else {
                // If the label triple does not exist, then create a new global view.
                auto temp_it = label_view_mapping_.emplace(label_triple, std::make_pair(std::vector<Edge>(), valid_view_count_));
                view_it = temp_it.first;
                view_it->second.first.push_back(e);
                edge_view_mapping_[e] = view_it->second.second;

                // Create a global view.
                auto &view = views_[valid_view_count_++];

                uint32_t u_candidates_count;
                auto u_candidates = data_graph->getVerticesByLabel(u_label, u_candidates_count);

                for (uint32_t j = 0; j < u_candidates_count; ++j) {
                    uint32_t v = u_candidates[j];
                    auto it = view.trie_.emplace(v, std::vector<uint32_t>());
                    auto &nbrs = it.first->second;

                    uint32_t v_nbrs_count;
                    auto v_nbrs = data_graph->getVertexNeighbors(v, v_nbrs_count);

                    for (uint32_t k = 0; k < v_nbrs_count; ++k) {
                        uint32_t vv = v_nbrs[k];
                        uint32_t vv_label = data_graph->getVertexLabel(vv);
                        uint32_t data_edge_label = data_graph->getEdgeLabelByLocalOffset(v, k);

                        if (vv_label == uu_label && query_edge_label == data_edge_label) {
                            std::pair<LabelID, LabelID> key = {data_edge_label, vv_label};
                            uint32_t nlf_id =  nlf_mapping_[key];
                            data_nlf_[v][nlf_id] += 1;
                            nbrs.push_back(vv);
                        }
                    }

                    view.cardinality_ += nbrs.size();
                }
            }
        }
    }

    // Initialize the global candidates for each query vertex.
    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            uint32_t u_label = query_graph->getVertexLabel(u);
            uint32_t u_candidates_count;
            auto u_candidates = data_graph->getVerticesByLabel(u_label, u_candidates_count);

            for (uint32_t i = 0; i < u_candidates_count; ++i) {
                uint32_t v = u_candidates[i];
                if (nlf_check(u, v)) {
                    uint32_t encoded_id = encoded_candidates_[u].size();
                    candidates_[u][v] = encoded_id;
                    encoded_candidates_[u].push_back(v);
                }
            }
        }
    }

    // Initialize the nlf views for each non leaf edge.
    for (uint32_t u = 0; u < query_graph->getVerticesCount(); ++u) {
        if (query_graph->getVertexDegree(u) > 1) {
            uint32_t u_nbrs_count;
            auto u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);

            for (uint32_t i = 0; i < u_nbrs_count; ++i) {
                uint32_t uu = u_nbrs[i];
                if (query_graph->getVertexDegree(uu) > 1) {
                    Edge query_edge = {u, uu};
                    edge_nlf_view_mapping_[query_edge] = valid_nlf_view_count_;
                    auto &view = nlf_views_[valid_nlf_view_count_++];
                    auto gv = get_view(query_edge);
                    auto &encoded_u_candidates = encoded_candidates_[u];
                    auto &uu_candidates = candidates_[uu];

                    for (auto v : encoded_u_candidates) {
                        view.adj_list_.emplace_back(std::vector<uint32_t>());
                        auto& nbrs = view.adj_list_.back();

                        uint32_t v_nbrs_count;
                        auto v_nbrs = gv->get_neighbor(v, v_nbrs_count);

                        for (uint32_t j = 0; j < v_nbrs_count; ++j) {
                            uint32_t vv = v_nbrs[j];
                            if (uu_candidates.contains(vv)) {
                                uint32_t encoded_vv = uu_candidates[vv];
                                nbrs.push_back(encoded_vv);
                            }
                        }

                        std::sort(nbrs.begin(), nbrs.end());
                        view.cardinality_ += nbrs.size();
                    }
                }
            }
        }
    }
}

MappedViews *GlobalViewManager::get_mapped_views(LabelTriple label_triple) {
    auto it = label_view_mapping_.find(label_triple);
    if (it != label_view_mapping_.end()) {
        return &it->second;
    }
    return nullptr;
}

GlobalEdgeView *GlobalViewManager::get_view(uint32_t view_id) {
    assert(view_id < valid_view_count_);
    return &views_[view_id];
}


GlobalEdgeView *GlobalViewManager::get_view(Edge edge) {
    auto it = edge_view_mapping_.find(edge);
    if (it != edge_view_mapping_.end())
        return &views_[it->second];
    return nullptr;
}

EncodedEdgeView *GlobalViewManager::get_nlf_view(Edge edge) {
    auto it = edge_nlf_view_mapping_.find(edge);
    if (it != edge_nlf_view_mapping_.end())
        return &nlf_views_[it->second];
    return nullptr;
}

uint32_t GlobalViewManager::get_view_id(Edge edge) {
    auto it = edge_view_mapping_.find(edge);
    if (it != edge_view_mapping_.end()) {
        return it->second;
    }
    return valid_view_count_;
}

uint32_t GlobalViewManager::get_encoded_candidate_set_size(uint32_t u) {
    return encoded_candidates_[u].size();
}

void GlobalViewManager::update_view(char op, uint32_t view_id, Edge e, LabelTriple label_triple) {
    auto& view = views_[view_id];
    std::pair<LabelID, LabelID> key = {label_triple.edge_label_, label_triple.dst_label_};
    uint32_t id = nlf_mapping_[key];

    if (op == '+') {
        if (view.insert(e.first, e.second)) {
            data_nlf_[e.first][id] += 1;
            // Update candidates
            auto mapped_edges = get_mapped_views(label_triple);
            for (auto& query_edge : mapped_edges->first) {
                if (query_graph_->getVertexDegree(query_edge.first) > 1) {
                    if (nlf_check(query_edge.first, e.first)) {
                        if (!candidates_[query_edge.first].contains(e.first)) {
                            uint32_t encoded_id = encoded_candidates_[query_edge.first].size();
                            candidates_[query_edge.first][e.first] = encoded_id;
                            encoded_candidates_[query_edge.first].push_back(e.first);
                            updated_candidates_.emplace_back(query_edge.first, e.first);
                        }
                    }
                }
            }
        }
    }
    else {
        if (view.remove(e.first, e.second)) {
            data_nlf_[e.first][id] -= 1;
            // Update candidates
            auto mapped_edges = get_mapped_views(label_triple);
            for (auto& query_edge : mapped_edges->first) {
                if (query_graph_->getVertexDegree(query_edge.first) > 1) {
                    if (!nlf_check(query_edge.first, e.first)) {
                        if (candidates_[query_edge.first].contains(e.first)) {
                            // candidates_[query_edge.first].erase(e.first);
                            // TODO: record the hole in the recorded list.
                            updated_candidates_.emplace_back(query_edge.first, e.first);
                        }
                    }
                }
            }
        }
    }
}

void GlobalViewManager::update_nlf_view(char op, Edge e, LabelTriple label_triple) {
    // Update candidates. Update each reverse view and generate the candidate set.
    for (auto kv : updated_candidates_) {
        uint32_t u = kv.first;
        uint32_t v = kv.second;

        if (!candidates_[u].contains(v))
            continue;

        uint32_t u_nbrs_count;
        auto u_nbrs = query_graph_->getVertexNeighbors(u, u_nbrs_count);

        for (uint32_t i = 0; i < u_nbrs_count; ++i) {
            uint32_t uu = u_nbrs[i];
            if (query_graph_->getVertexDegree(uu) > 1) {
                Edge query_edge = {u, uu};
                Edge reverse_edge = {uu, u};
                auto reverse_nlf_view = get_nlf_view(reverse_edge);
                auto nlf_view = get_nlf_view(query_edge);

                if (op == '+') {
                    auto gv = get_view(query_edge);
                    uint32_t v_nbrs_count;
                    auto v_nbrs = gv->get_neighbor(v, v_nbrs_count);
                    uint32_t encoded_v = candidates_[u][v];
                    nlf_view->adj_list_.emplace_back(std::vector<uint32_t>());
                    auto &nbrs = nlf_view->adj_list_.back();

                    for (uint32_t j = 0; j < v_nbrs_count; ++j) {
                        uint32_t vv = v_nbrs[j];
                        if (candidate_check(uu, vv)) {
                            uint32_t encoded_vv = candidates_[uu][vv];

                            if (reverse_nlf_view->contains(encoded_vv))
                                reverse_nlf_view->insert(encoded_vv, encoded_v);

                            // Update the neighbors of v.
                            nbrs.push_back(encoded_vv);
                        }
                    }

                    std::sort(nbrs.begin(), nbrs.end());
                }
                else {
                    uint32_t encoded_v = candidates_[u][v];
                    auto &nbrs = nlf_view->adj_list_[encoded_v];

                    for (auto encoded_vv : nbrs) {
                        if (reverse_nlf_view->contains(encoded_vv))
                            reverse_nlf_view->remove(encoded_vv, encoded_v);
                    }

                    nbrs.clear();
                }
            }
        }

        if (op == '-') {
            candidates_[u].erase(v);
        }
    }

    auto mapped_views = get_mapped_views(label_triple);
    for (auto query_edge : mapped_views->first) {
        uint32_t u0 = query_edge.first;
        uint32_t u1 = query_edge.second;
        if (query_graph_->getVertexDegree(u0) > 1 && query_graph_->getVertexDegree(u1) > 1) {
            if (candidate_check(u0, e.first) && candidate_check(u1, e.second)) {
                auto nlf_view = get_nlf_view(query_edge);
                Edge reverse_query_edge = {u1, u0};
                auto reverse_nlf_view = get_nlf_view(reverse_query_edge);

                uint32_t encoded_v = candidates_[query_edge.first][e.first];
                uint32_t encoded_vv = candidates_[query_edge.second][e.second];

                if (op == '+') {
                    nlf_view->insert(encoded_v, encoded_vv);
                    reverse_nlf_view->insert(encoded_vv, encoded_v);
                }
                else {
                    nlf_view->remove(encoded_v, encoded_vv);
                    reverse_nlf_view->remove(encoded_vv, encoded_v);
                }
            }
        }
    }

    updated_candidates_.clear();
}

void GlobalViewManager::calculate_memory_usage() {

}

void GlobalViewManager::print_view_info() {
    printf("Global View Info:\n");
    for (auto& kv : edge_view_mapping_) {
        printf("e(%u, %u): vertex num %u, edge num %u\n", kv.first.first, kv.first.second,
               views_[kv.second].get_key_num(), views_[kv.second].get_edge_num());
    }

    printf("NLF View Info:\n");
    for (auto& kv : edge_nlf_view_mapping_) {
        printf("e(%u, %u): vertex num %u, edge num %u\n", kv.first.first, kv.first.second,
               nlf_views_[kv.second].get_key_num(), nlf_views_[kv.second].get_edge_num());
    }

    for (auto&kv : label_view_mapping_) {
        printf("(%u, %u, %u) ", kv.first.src_label_, kv.first.edge_label_, kv.first.dst_label_);
    }
    printf("\n");
    for (auto& candidates : candidates_)
        printf("%zu, ", candidates.size());

    printf("\n");
}