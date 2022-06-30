//
// Created by ssunah on 6/22/18.
//

#ifndef SUBGRAPHMATCHING_GRAPH_H
#define SUBGRAPHMATCHING_GRAPH_H

#include <unordered_map>
#include <iostream>
#include <vector>
#include "utility/sparsepp/spp.h"
#include "configuration/types.h"
#include "configuration/config.h"

/**
 * A graph is stored as the CSR format.
 */
using spp::sparse_hash_map;
class Graph {
private:
    bool enable_label_offset_;

    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;

    ui* offsets_;
    VertexID * neighbors_;
    LabelID* vertex_labels_;
    ui* reverse_index_offsets_;
    ui* reverse_index_;

    LabelID *edge_labels_;

    int* core_table_;
    ui core_length_;

    std::unordered_map<LabelID, ui> labels_frequency_;
    sparse_hash_map<uint64_t, std::vector<edge>* >* edge_index_;

#if OPTIMIZED_LABELED_GRAPH == 1
    ui* labels_offsets_;
    sparse_hash_map<uint32_t, uint32_t>* nlf_;
#endif

private:
    void BuildReverseIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    void BuildNLF();
    void BuildLabelOffset();
#endif

public:
    bool is_edge_labeled;
public:
    Graph(const bool enable_label_offset) {
        enable_label_offset_ = enable_label_offset;

        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;
        core_length_ = 0;

        is_edge_labeled = false;
        offsets_ = NULL;
        neighbors_ = NULL;
        vertex_labels_ = NULL;
        reverse_index_offsets_ = NULL;
        reverse_index_ = NULL;
        core_table_ = NULL;
        labels_frequency_.clear();
        edge_index_ = NULL;
        edge_labels_ = NULL;
#if OPTIMIZED_LABELED_GRAPH == 1
        labels_offsets_ = NULL;
        nlf_ = NULL;
#endif
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] vertex_labels_;
        delete[] reverse_index_offsets_;
        delete[] reverse_index_;
        delete[] core_table_;
        delete[] edge_labels_;
        delete edge_index_;
#if OPTIMIZED_LABELED_GRAPH == 1
        delete[] labels_offsets_;
        delete[] nlf_;
#endif
    }

public:
    void loadGraphFromMemory(std::vector<std::pair<uint32_t, uint32_t>> &vertex_list,
                             std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &edge_list);

    void loadGraphFromFile(const std::string& file_path);

    void loadGraphFromFileWithoutMeta(const std::string &file_path);

    void loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                     const std::string &label_path, const std::string &edge_label_path);

    void loadGraphFromFileCompressed(const std::string &file_path);

    void storeCompressedGraph(const std::string &degree_path, const std::string &edge_path,
                              const std::string &label_path, const std::string &edge_label_file_path);

    void storeCompressedGraph(const std::string &file_path);

    void printGraphMetaData();
public:
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }

    const ui getCoreValue(const VertexID id) const {
        return core_table_[id];
    }

    const ui get2CoreSize() const {
        return core_length_;
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return vertex_labels_[id];
    }

    LabelID getEdgeLabel(ui edge_id) const {
        return edge_labels_[edge_id];
    }

    LabelID getEdgeLabelByLocalOffset(ui id, ui offset) const {
        return getEdgeLabel(offset + offsets_[id]);
    }

    uint32_t getEdgeLabelByVertex(uint32_t u, uint32_t v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int mid = 0;
        int begin = 0;
        int end = count - 1;
        bool find = false;
        while (begin <= end) {
            mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                find = true;
                break;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return find ? getEdgeLabelByLocalOffset(v, mid) : std::numeric_limits<uint32_t>::max();
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    const sparse_hash_map<uint64_t, std::vector<edge>*>* getEdgeIndex() const {
        return edge_index_;
    }

    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets_[id + 1] - reverse_index_offsets_[id];
        return reverse_index_ + reverse_index_offsets_[id];
    }

#if OPTIMIZED_LABELED_GRAPH == 1
    const ui * getNeighborsByLabel(const VertexID id, const LabelID label, ui& count) const {
        ui offset = id * labels_count_ + label;
        count = labels_offsets_[offset + 1] - labels_offsets_[offset];
        return neighbors_ + labels_offsets_[offset];
    }

    const sparse_hash_map<uint32_t, uint32_t>* getVertexNLF(const VertexID id) const {
        return nlf_ + id;
    }

    bool checkEdgeExistence(const VertexID u, const VertexID v, const LabelID u_label) const {
        ui count = 0;
        const VertexID* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
#endif

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void buildCoreTable();

    void buildEdgeIndex();
};


#endif //SUBGRAPHMATCHING_GRAPH_H
