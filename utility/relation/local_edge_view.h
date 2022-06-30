//
// Created by sunsx on 28/05/21.
//

#ifndef RAPIDMATCH_LOCALEDGEVIEW_CPP
#define RAPIDMATCH_LOCALEDGEVIEW_CPP

#include <vector>
#include <algorithm>
#include "../sparsepp/spp.h"
using spp::sparse_hash_map;

class LocalEdgeView {
public:
    // Key: vertex id; Value: first is the begin position and the second is the end position.
    std::vector<std::pair<uint32_t, uint32_t>> trie_;

    // The data array.
    uint32_t* data_;

    // The number of edges in the view.
    uint32_t cardinality_;
public:
    LocalEdgeView() {
        trie_.reserve(1024);
        cardinality_ = 0;
    }

    ~LocalEdgeView() {
    }

    uint32_t get_neighbor_num(uint32_t key) {
        return trie_[key].second - trie_[key].first;
    }

    bool contains(uint32_t key) {
        return key < trie_.size();
    }

    uint32_t* get_neighbors(uint32_t key, uint32_t& count) {
        count = trie_[key].second - trie_[key].first;
        return data_ + trie_[key].first;
    }

    uint32_t get_key_num() {
        return static_cast<uint32_t>(trie_.size());
    }
    uint32_t get_edge_num() {
        return cardinality_;
    }

    void clear() {
        trie_.clear();
        data_ = nullptr;
        cardinality_ = 0;
    }

    uint64_t memory_cost() {
        if (get_key_num() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        double hash_table_memory_cost = static_cast<double>((get_key_num() * 3 * per_element_size));
        memory_cost += get_edge_num() * per_element_size + static_cast<uint64_t>(hash_table_memory_cost);
        return memory_cost;
    }
};

#endif //RAPIDMATCH_LOCALEDGEVIEW_CPP
