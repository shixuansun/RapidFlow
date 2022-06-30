//
// Created by sunsx on 28/05/21.
//

#ifndef RAPIDMATCH_STREAMING_TYPE_H
#define RAPIDMATCH_STREAMING_TYPE_H

#include <cstdint>
#include <vector>

#include "sparsepp/spp.h"


/**
 * begin label, edge label, end label;
 */
typedef struct LabelTriple {
    uint32_t src_label_;
    uint32_t edge_label_;
    uint32_t dst_label_;
    bool operator==(const LabelTriple& l) const {
        return l.src_label_ == src_label_ && l.edge_label_ == edge_label_ && l.dst_label_ == dst_label_;
    }
} LabelTriple;

/**
 * begin vertex, end vertex
 */
typedef std::pair<uint32_t, uint32_t> Edge;

/**
 * Edges and the corresponding view.
 */
typedef std::pair<std::vector<Edge>, uint32_t> MappedViews;

/**
 * Orders for each edge
 */

enum RelationEdgeType {
    REGULAR = 0,
    FULL_CONNECTION = 1
};

typedef struct OrdersPerEdge {
    std::vector<uint32_t> indexing_order_;
    uint32_t triangle_end_index_;
    uint32_t adjacent_end_index_;
    std::vector<uint32_t> indexing_order_bn_offset_;
    std::vector<uint32_t> indexing_order_bn_;
    std::vector<uint32_t> matching_order_;
    std::vector<uint32_t> matching_order_bn_offset_;
    std::vector<uint32_t> matching_order_bn_;
    std::vector<uint32_t> matching_order_view_mappings_;
    std::vector<RelationEdgeType> matching_order_edge_type_;
} OrdersPerEdge;

namespace std
{
    template<>
    struct hash<LabelTriple>{
        std::size_t operator()(LabelTriple const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.src_label_);
            spp::hash_combine(seed, l.dst_label_);
            spp::hash_combine(seed, l.edge_label_);
            return seed;
        }
    };

    template<>
    struct hash<Edge>{
        std::size_t operator()(Edge const &l) const{
            std::size_t seed = 0;
            spp::hash_combine(seed, l.first);
            spp::hash_combine(seed, l.second);
            return seed;
        }
    };
}

typedef struct Update {
    uint64_t id_;
    char op_;
    Edge edge_;
    LabelTriple labels_;
} Update;


typedef struct VertexOrderingPriority {
    uint32_t bn_count_;
    uint32_t core_value_;
    uint32_t degree_;
    uint32_t vertex_id_;
    bool operator <(const VertexOrderingPriority& rhs) const {
        if (bn_count_ != rhs.bn_count_)
            return bn_count_ < rhs.bn_count_;

        if (core_value_ != rhs.core_value_)
            return core_value_ < rhs.core_value_;

        if (degree_ != rhs.degree_)
            return degree_ < rhs.degree_;

        return vertex_id_ < rhs.vertex_id_;
    }
} VertexOrderingPriority;

#endif //RAPIDMATCH_STREAMING_TYPE_H
