#ifndef SUBGRAPHMATCHING_GRAPHOPERATIONS_H
#define SUBGRAPHMATCHING_GRAPHOPERATIONS_H

#include "graph/graph.h"
class GraphOperations {
public:
    static void getKCore(const Graph *graph, int *core_table);
    static void compute_degeneracy_order(const Graph* graph, uint32_t* degeneracy_order);
    static void compute_automorphism(const Graph* graph, std::vector<std::vector<uint32_t>>& embeddings);
private:
    static void dfs(TreeNode* tree, VertexID cur_vertex, VertexID* dfs_order, ui& count);
};


#endif //SUBGRAPHMATCHING_GRAPHOPERATIONS_H
