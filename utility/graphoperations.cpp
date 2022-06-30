#include "graphoperations.h"
#include <memory.h>
#include <queue>

void GraphOperations::getKCore(const Graph *graph, int *core_table) {
    int vertices_count = graph->getVerticesCount();
    int max_degree = graph->getGraphMaxDegree();

    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Degree from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        ui count;
        const VertexID * neighbors = graph->getVertexNeighbors(v, count);

        for(int j = 0; j < static_cast<int>(count); ++j) {
            int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }
    }

    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}


void GraphOperations::dfs(TreeNode *tree, VertexID cur_vertex, VertexID *dfs_order, ui &count) {
    dfs_order[count++] = cur_vertex;

    for (ui i = 0; i < tree[cur_vertex].children_count_; ++i) {
        dfs(tree, tree[cur_vertex].children_[i], dfs_order, count);
    }
}

void GraphOperations::compute_degeneracy_order(const Graph *graph, uint32_t *degeneracy_order) {
    int vertices_count = graph->getVerticesCount();
    int max_degree = graph->getGraphMaxDegree();

    int* core_table = new int[vertices_count];        // core values.
    int* vertices = new int[vertices_count];          // Vertices sorted by degree.
    int* position = new int[vertices_count];          // The position of vertices in vertices array.
    int* degree_bin = new int[max_degree + 1];      // Degree from 0 to max_degree.
    int* offset = new int[max_degree + 1];          // The offset in vertices array according to degree.

    std::fill(degree_bin, degree_bin + (max_degree + 1), 0);

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        core_table[i] = degree;
        degree_bin[degree] += 1;
    }

    int start = 0;
    for (int i = 0; i < max_degree + 1; ++i) {
        offset[i] = start;
        start += degree_bin[i];
    }

    for (int i = 0; i < vertices_count; ++i) {
        int degree = graph->getVertexDegree(i);
        position[i] = offset[degree];
        vertices[position[i]] = i;
        offset[degree] += 1;
    }

    for (int i = max_degree; i > 0; --i) {
        offset[i] = offset[i - 1];
    }
    offset[0] = 0;

    for (int i = 0; i < vertices_count; ++i) {
        int v = vertices[i];

        ui count;
        const VertexID * neighbors = graph->getVertexNeighbors(v, count);

        for(int j = 0; j < static_cast<int>(count); ++j) {
            int u = neighbors[j];

            if (core_table[u] > core_table[v]) {

                // Get the position and vertex which is with the same degree
                // and at the start position of vertices array.
                int cur_degree_u = core_table[u];
                int position_u = position[u];
                int position_w = offset[cur_degree_u];
                int w = vertices[position_w];

                if (u != w) {
                    // Swap u and w.
                    position[u] = position_w;
                    position[w] = position_u;
                    vertices[position_u] = w;
                    vertices[position_w] = u;
                }

                offset[cur_degree_u] += 1;
                core_table[u] -= 1;
            }
        }

        degeneracy_order[i] = v;
    }

    delete[] core_table;
    delete[] vertices;
    delete[] position;
    delete[] degree_bin;
    delete[] offset;
}

void GraphOperations::compute_automorphism(const Graph *graph, std::vector<std::vector<uint32_t>> &embeddings) {
    // Note that this method is working on small graphs (tens of vertices) only.
    // Initialize resource.
    uint32_t n = graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    std::vector<uint32_t> idx(n);
    std::vector<uint32_t> mapping(n);
    std::vector<std::vector<uint32_t>> local_candidates(n);
    std::vector<std::vector<uint32_t>> global_candidates(n);
    std::vector<std::vector<uint32_t>> backward_neighbors(n);

    // Initialize global candidates.
    for (uint32_t u = 0; u < n; ++u) {
        uint32_t u_label = graph->getVertexLabel(u);
        uint32_t u_degree = graph->getVertexDegree(u);

        for (uint32_t v = 0; v < n; ++v) {
            uint32_t v_label = graph->getVertexLabel(v);
            uint32_t v_degree = graph->getVertexDegree(v);

            if (v_label == u_label && v_degree >= u_degree)
                global_candidates[u].push_back(v);
        }
    }

    // Generate a matching order.
    std::vector<uint32_t> matching_order;
    uint32_t selected_vertex = 0;
    uint32_t selected_vertex_selectivity = global_candidates[selected_vertex].size();
    for (uint32_t u = 1; u < n; ++u) {
        if (global_candidates[u].size() < selected_vertex_selectivity){
            selected_vertex = u;
            selected_vertex_selectivity = global_candidates[u].size();
        }
    }

    matching_order.push_back(selected_vertex);
    visited[selected_vertex] = true;

    for (uint32_t i = 1; i < n; ++i) {
        selected_vertex_selectivity = n + 1;
        for (uint32_t u = 0; u < n; ++u) {
            if (!visited[u]) {
                bool is_feasible = false;

                uint32_t u_nbr_count;
                auto u_nbr = graph->getVertexNeighbors(u, u_nbr_count);
                for (uint32_t j = 0; j < u_nbr_count; ++j) {
                    uint32_t uu = u_nbr[j];

                    if (visited[uu]) {
                        is_feasible = true;
                        break;
                    }
                }

                if (is_feasible && global_candidates[u].size() < selected_vertex_selectivity) {
                    selected_vertex = u;
                    selected_vertex_selectivity = global_candidates[u].size();
                }
            }
        }
        matching_order.push_back(selected_vertex);
        visited[selected_vertex] = true;
    }

    std::fill(visited.begin(), visited.end(), false);

    // Set backward neighbors to compute local candidates.
    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = matching_order[i];
        for (uint32_t j = 0; j < i; ++j) {
            uint32_t uu = matching_order[j];

            if (graph->checkEdgeExistence(uu, u)) {
                backward_neighbors[u].push_back(uu);
            }
        }
    }

    // Recursive search along the matching order.
    int cur_level = 0;
    local_candidates[cur_level] = global_candidates[matching_order[0]];

    while (true) {
        while (idx[cur_level] < local_candidates[cur_level].size()) {
            uint32_t u = matching_order[cur_level];
            uint32_t v = local_candidates[cur_level][idx[cur_level]];
            idx[cur_level] += 1;

            if (cur_level == n - 1) {
                // Find an embedding
                mapping[u] = v;
                embeddings.push_back(mapping);
            } else {
                mapping[u] = v;
                visited[v] = true;
                cur_level += 1;
                idx[cur_level] = 0;

                {
                    // Compute local candidates.
                    u = matching_order[cur_level];
                    for (auto temp_v: global_candidates[u]) {
                        if (!visited[temp_v]) {
                            bool is_feasible = true;
                            for (auto uu: backward_neighbors[u]) {
                                uint32_t edge_label = graph->getEdgeLabelByVertex(u, uu);

                                uint32_t temp_vv = mapping[uu];
                                // TODO: check edge label
                                if (!graph->checkEdgeExistence(temp_v, temp_vv)
                                    || edge_label != graph->getEdgeLabelByVertex(temp_v, temp_vv)) {
                                    is_feasible = false;
                                    continue;
                                }
                            }
                            if (is_feasible)
                                local_candidates[cur_level].push_back(temp_v);
                        }
                    }
                }
            }
        }

        local_candidates[cur_level].clear();
        cur_level -= 1;
        if (cur_level < 0) {
            break;
        }
        visited[mapping[matching_order[cur_level]]] = false;
    }
}

