//
// Created by sunsx on 28/05/21.
//

#ifndef RAPIDMATCH_ORDER_MANAGER_H
#define RAPIDMATCH_ORDER_MANAGER_H

#include "graph/graph.h"
#include "streaming_type.h"

class OrderManager {
private:
    std::vector<std::vector<uint32_t>> automorphisms_;
    std::vector<OrdersPerEdge> orders_;
    spp::sparse_hash_map<Edge, uint32_t> edge_orders_mapping_;
    spp::sparse_hash_map<LabelTriple, std::vector<Edge>> label_edge_mapping_;
    spp::sparse_hash_map<LabelTriple, std::vector<uint32_t>> label_automorphism_mapping_;
    std::vector<std::vector<Edge>> automorphism_edges_;

private:
    // New framework.
    void detect_automorphism_edges(const Graph *query_graph);
    void create_indexing_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    void create_matching_order_for_reduced_graph(const Graph *query_graph, Edge edge, OrdersPerEdge &order);
    void generate_matching_order_with_RI(const Graph* graph, std::vector<uint32_t>& matching_order);

public:
    OrderManager() {}
    ~OrderManager() { release(); }

    void initialize(const Graph *query_graph);
    void release();

    OrdersPerEdge* get_orders(Edge edge);
    std::vector<Edge>* get_mapped_edges(LabelTriple label_triple);
    std::vector<uint32_t>* get_mapped_automorphism(LabelTriple labelTriple);
    std::tuple<Edge, OrdersPerEdge*, uint32_t> get_automorphism_meta(uint32_t id);

    void create_orders(const Graph *query_graph);
    void print_info();
};


#endif //RAPIDMATCH_ORDER_MANAGER_H
