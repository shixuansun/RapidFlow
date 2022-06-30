//
// Created by sunsh on 11/21/2021.
//

#ifndef RAPIDMATCH_ENCODED_EDGE_VIEW_H
#define RAPIDMATCH_ENCODED_EDGE_VIEW_H

#include <vector>

class EncodedEdgeView {
public:
    // Key: a vertex; Value: neighbors of the vertex in key.
    std::vector<std::vector<uint32_t>> adj_list_;
    uint32_t cardinality_;

public:
    EncodedEdgeView() { cardinality_ = 0; }

    ~EncodedEdgeView() {}

    bool insert(uint32_t u, uint32_t v) {
        auto& neighbors = adj_list_[u];
        auto lb = std::lower_bound(neighbors.begin(), neighbors.end(), v);
        if (lb == neighbors.end() || *lb != v) {
            neighbors.insert(lb, v);
            cardinality_ += 1;
            return true;
        }
        return false;
    }

    bool remove(uint32_t u, uint32_t v) {
        auto& neighbors = adj_list_[u];
        auto lb = std::lower_bound(neighbors.begin(), neighbors.end(), v);
        if (lb != neighbors.end() && *lb == v) {
            neighbors.erase(lb);
            cardinality_ -= 1;
            return true;
        }
        return false;
    }

    uint32_t get_neighbor_num(uint32_t key) {
        return adj_list_[key].size();
    }

    uint32_t* get_neighbor(uint32_t key, uint32_t& count) {
        count = adj_list_[key].size();
        return adj_list_[key].data();
    }

    bool contains(uint32_t key) {
        return key < adj_list_.size();
    }

    uint32_t get_key_num() {
        return static_cast<uint32_t>(adj_list_.size());
    }

    uint32_t get_edge_num() {
        return cardinality_;
    }

    uint64_t memory_cost() {
        if (get_key_num() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_key_num() * (per_element_size + sizeof(std::vector<uint32_t>)) + get_edge_num() * per_element_size;
        return memory_cost;
    }
};

#endif //RAPIDMATCH_ENCODED_EDGE_VIEW_H
