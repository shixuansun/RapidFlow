//
// Created by sunsx on 5/28/21.
//

#ifndef INTERNALRAPIDMATCH_VIEW_H
#define INTERNALRAPIDMATCH_VIEW_H

#include <vector>
#include <algorithm>
#include "../sparsepp/spp.h"
using spp::sparse_hash_map;

class GlobalEdgeView {
public:
    // Key: a vertex; Value: neighbors of the vertex in key.
    sparse_hash_map<uint32_t, std::vector<uint32_t>> trie_;
    uint32_t cardinality_;

public:
    GlobalEdgeView() { cardinality_ = 0; }

    ~GlobalEdgeView() {}

    bool insert(uint32_t u, uint32_t v) {
        auto it = trie_.find(u);
        if (it != trie_.end()) {
            auto lb = std::lower_bound(it->second.begin(), it->second.end(), v);
            if (lb == it->second.end() || *lb != v) {
                it->second.insert(lb, v);
                cardinality_ += 1;
                return true;
            }
            return false;
        }
        else {
            auto temp_it = trie_.emplace(u, std::vector<uint32_t>());
            temp_it.first->second.push_back(v);
            cardinality_ += 1;
            return true;
        }
    }

    bool remove(uint32_t u, uint32_t v) {
        auto iter = trie_.find(u);
        if (iter != trie_.end()) {
            auto lb = std::lower_bound(iter->second.begin(), iter->second.end(), v);
            if (lb != iter->second.end() && *lb == v) {
                iter->second.erase(lb);
                cardinality_ -= 1;
                return true;
            }
        }

        return false;
    }

    uint32_t get_neighbor_num(uint32_t key) {
        auto iter = trie_.find(key);
        if (iter != trie_.end())
            return iter->second.size();
        return 0;
    }

    uint32_t* get_neighbor(uint32_t key, uint32_t& count) {
        count = 0;
        auto iter = trie_.find(key);
        if (iter != trie_.end()) {
            count = iter->second.size();
            return iter->second.data();
        }
        return nullptr;
    }

    bool contains(uint32_t key) {
        return trie_.contains(key);
    }

    uint32_t get_key_num() {
        return static_cast<uint32_t>(trie_.size());
    }

    uint32_t get_edge_num() {
        return cardinality_;
    }

    uint64_t memory_cost() {
        if (get_key_num() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_key_num() * (per_element_size + sizeof(std::vector<uint32_t>)) * trie_.load_factor() + get_edge_num() * per_element_size;
        return memory_cost;
    }
};


#endif //INTERNALRAPIDMATCH_VIEW_H
