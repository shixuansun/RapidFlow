//
// Created by ssunah on 6/22/18.
//

#include "graph.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <utility/graphoperations.h>

const char* degree_file_name = "b_degree.bin";
const char* adj_file_name = "b_adj.bin";
const char* vertex_label_file_name = "b_vertex_label.bin";
const char* edge_label_file_name = "b_edge_label.bin";

void Graph::BuildReverseIndex() {
    reverse_index_ = new ui[vertices_count_];
    reverse_index_offsets_= new ui[labels_count_ + 1];
    reverse_index_offsets_[0] = 0;

    ui total = 0;
    for (ui i = 0; i < labels_count_; ++i) {
        reverse_index_offsets_[i + 1] = total;
        total += labels_frequency_[i];
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID label = vertex_labels_[i];
        reverse_index_[reverse_index_offsets_[label + 1]++] = i;
    }
}

#if OPTIMIZED_LABELED_GRAPH == 1
void Graph::BuildNLF() {
    nlf_ = new sparse_hash_map<uint32_t, uint32_t>[vertices_count_];
    for (ui i = 0; i < vertices_count_; ++i) {
        ui count;
        const VertexID * neighbors = getVertexNeighbors(i, count);

        for (ui j = 0; j < count; ++j) {
            VertexID u = neighbors[j];
            LabelID label = getVertexLabel(u);
            nlf_[i][label] += 1;
        }
    }
}

void Graph::BuildLabelOffset() {
    size_t labels_offset_size = (size_t)vertices_count_ * labels_count_ + 1;
    labels_offsets_ = new ui[labels_offset_size];
    std::fill(labels_offsets_, labels_offsets_ + labels_offset_size, 0);

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1],
            [this](const VertexID u, const VertexID v) -> bool {
                return vertex_labels_[u] == vertex_labels_[v] ? u < v : vertex_labels_[u] < vertex_labels_[v];
            });
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        LabelID previous_label = 0;
        LabelID current_label = 0;

        labels_offset_size = i * labels_count_;
        labels_offsets_[labels_offset_size] = offsets_[i];

        for (ui j = offsets_[i]; j < offsets_[i + 1]; ++j) {
            current_label = vertex_labels_[neighbors_[j]];

            if (current_label != previous_label) {
                for (ui k = previous_label + 1; k <= current_label; ++k) {
                    labels_offsets_[labels_offset_size + k] = j;
                }
                previous_label = current_label;
            }
        }

        for (ui l = current_label + 1; l <= labels_count_; ++l) {
            labels_offsets_[labels_offset_size + l] = offsets_[i + 1];
        }
    }
}

#endif

void Graph::loadGraphFromMemory(std::vector<std::pair<uint32_t, uint32_t>> &vertex_list,
                                std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> &edge_list) {
    vertices_count_ = vertex_list.size();
    edges_count_ = edge_list.size();

    offsets_ = new uint32_t[vertices_count_ + 1];
    std::fill(offsets_, offsets_ + vertices_count_ + 1, 0);

    neighbors_ = new VertexID[edges_count_ * 2];
    vertex_labels_ = new LabelID[vertices_count_];

    if (is_edge_labeled)
        edge_labels_ = new LabelID[edges_count_ * 2];

    labels_count_ = 0;
    max_degree_ = 0;

    uint32_t max_label_id = 0;
    for (auto& vertex : vertex_list) {
        uint32_t id = vertex.first;
        uint32_t label = vertex.second;
        vertex_labels_[id] = label;

        if (labels_frequency_.find(label) == labels_frequency_.end()) {
            labels_frequency_[label] = 0;
            if (label > max_label_id)
                max_label_id = label;
        }
        labels_frequency_[label] += 1;
    }

    for (auto& edge : edge_list) {
        uint32_t begin = std::get<0>(edge);
        uint32_t end = std::get<1>(edge);
        offsets_[begin] += 1;
        offsets_[end] += 1;
    }
    uint32_t temp = offsets_[0];
    offsets_[0] = 0;
    for (uint32_t i = 1; i <= vertices_count_; ++i) {
        uint32_t sum = temp + offsets_[i - 1];
        temp = offsets_[i];
        offsets_[i] = sum;
    }

    for (uint32_t i = 0; i < vertices_count_; ++i) {
        uint32_t degree = offsets_[i + 1] - offsets_[i];
        if (degree > max_degree_)
            max_degree_ = degree;
    }

    std::vector<uint32_t> local_offset(vertices_count_, 0);
    for (auto& edge : edge_list) {
        uint32_t begin = std::get<0>(edge);
        uint32_t end = std::get<1>(edge);
        uint32_t label = std::get<2>(edge);

        uint32_t offset = offsets_[begin] + local_offset[begin];
        neighbors_[offset] = end;
        edge_labels_[offset] = label;

        offset = offsets_[end] + local_offset[end];
        neighbors_[offset] = begin;
        edge_labels_[offset] = label;

        local_offset[begin] += 1;
        local_offset[end] += 1;
    }

    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    BuildReverseIndex();
    buildEdgeIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        BuildLabelOffset();
    }
#endif
}

void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count_ >> edges_count_;
    offsets_ = new ui[vertices_count_ +  1];
    offsets_[0] = 0;

    neighbors_ = new VertexID[edges_count_ * 2];
    vertex_labels_ = new LabelID[vertices_count_];

    if (is_edge_labeled)
        edge_labels_ = new LabelID[edges_count_ * 2];

    labels_count_ = 0;
    max_degree_ = 0;

    LabelID max_label_id = 0;
    std::vector<ui> neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            vertex_labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + degree;

            if (degree > max_degree_) {
                max_degree_ = degree;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            LabelID edge_label;

            if (is_edge_labeled) {
                infile >> begin >> end >> edge_label;
                ui offset = offsets_[begin] + neighbors_offset[begin];
                neighbors_[offset] = end;
                edge_labels_[offset] = edge_label;

                offset = offsets_[end] + neighbors_offset[end];
                neighbors_[offset] = begin;
                edge_labels_[offset] = edge_label;
            }
            else {
                infile >> begin >> end;
                ui offset = offsets_[begin] + neighbors_offset[begin];
                neighbors_[offset] = end;

                offset = offsets_[end] + neighbors_offset[end];
                neighbors_[offset] = begin;
            }

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;
        }
    }

    infile.close();
    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }

    BuildReverseIndex();
    buildEdgeIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        // BuildLabelOffset();
    }
#endif
}

void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count_ << ", |E|: " << edges_count_ << ", |\u03A3|: " << labels_count_ << std::endl;
    std::cout << "Max Degree: " << max_degree_ << ", Max Label Frequency: " << max_label_frequency_ << std::endl;
}

void Graph::buildCoreTable() {
    core_table_ = new int[vertices_count_];
    GraphOperations::getKCore(this, core_table_);

    for (ui i = 0; i < vertices_count_; ++i) {
        if (core_table_[i] > 1) {
            core_length_ += 1;
        }
    }
}

void Graph::loadGraphFromFileCompressed(const std::string &degree_path, const std::string &edge_path,
                                        const std::string &label_path,
                                        const std::string &edge_label_path) {
    std::ifstream deg_file(degree_path, std::ios::binary);

    if (!deg_file.is_open()) {
        std::cerr << "Cannot open degree file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    auto start = std::chrono::high_resolution_clock::now();
    int int_size;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);
    deg_file.read(reinterpret_cast<char *>(&vertices_count_), 4);
    deg_file.read(reinterpret_cast<char *>(&edges_count_), 4);

    offsets_ = new ui[vertices_count_ + 1];
    ui* degrees = new unsigned int[vertices_count_];

    deg_file.read(reinterpret_cast<char *>(degrees), sizeof(int) * vertices_count_);


    deg_file.close();
    deg_file.clear();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Load degree file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    std::ifstream adj_file(edge_path, std::ios::binary);

    if (!adj_file.is_open()) {
        std::cerr << "Cannot open edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();
    size_t neighbors_count = (size_t)edges_count_ * 2;
    neighbors_ = new ui[neighbors_count];

    offsets_[0] = 0;
    for (ui i = 1; i <= vertices_count_; ++i) {
        offsets_[i] = offsets_[i - 1] + degrees[i - 1];
    }

    max_degree_ = 0;

    for (ui i = 0; i < vertices_count_; ++i) {
        if (degrees[i] > 0) {
            if (degrees[i] > max_degree_)
                max_degree_ = degrees[i];
            adj_file.read(reinterpret_cast<char *>(neighbors_ + offsets_[i]), degrees[i] * sizeof(int));
        }
    }

    adj_file.close();
    adj_file.clear();

    delete[] degrees;

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load adj file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    if (is_edge_labeled) {
        std::ifstream edge_label_file(edge_label_path, std::ios::binary);

        if (!edge_label_file.is_open()) {
            std::cerr << "Cannot open edge label file " << edge_label_path << " ." << std::endl;
            exit(-1);
        }

        start = std::chrono::high_resolution_clock::now();
        size_t edge_label_size = (size_t)edges_count_ * 2;
        edge_labels_ = new LabelID[edge_label_size];

        for (ui i = 0; i < vertices_count_; ++i) {
            ui degree = getVertexDegree(i);
            if (degree > 0) {
                edge_label_file.read(reinterpret_cast<char *>(edge_labels_ + offsets_[i]), degree * sizeof(int));
            }
        }

        edge_label_file.close();
        edge_label_file.clear();

        end = std::chrono::high_resolution_clock::now();
        std::cout << "Load edge label file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
    }

    std::ifstream label_file(label_path, std::ios::binary);
    if (!label_file.is_open())  {
        std::cerr << "Cannot open label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    start = std::chrono::high_resolution_clock::now();

    vertex_labels_ = new ui[vertices_count_];
    label_file.read(reinterpret_cast<char *>(vertex_labels_), sizeof(int) * vertices_count_);

    label_file.close();
    label_file.clear();

    ui max_label_id = 0;
    for (ui i = 0; i < vertices_count_; ++i) {
        ui label = vertex_labels_[i];

        if (labels_frequency_.find(label) == labels_frequency_.end()) {
            labels_frequency_[label] = 0;
            if (label > max_label_id)
                max_label_id = label;
        }

        labels_frequency_[label] += 1;
    }

    labels_count_ = (ui)labels_frequency_.size() > (max_label_id + 1) ? (ui)labels_frequency_.size() : max_label_id + 1;

    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Load label file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    BuildReverseIndex();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Build reverse index file time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    buildEdgeIndex();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Build edge index time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;

#if OPTIMIZED_LABELED_GRAPH == 1
    if (enable_label_offset_) {
        BuildNLF();
        BuildLabelOffset();
    }
#endif
}

void Graph::storeCompressedGraph(const std::string &degree_path, const std::string &edge_path,
                                 const std::string &label_path, const std::string &edge_label_file_path) {
    ui* degrees = new ui[vertices_count_];
    for (ui i = 0; i < vertices_count_; ++i) {
        degrees[i] = offsets_[i + 1] - offsets_[i];
    }

    std::ofstream deg_outputfile(degree_path, std::ios::binary);

    if (deg_outputfile.is_open()) {
        std::cout << "Open degree file " << degree_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot degree edge file " << degree_path << " ." << std::endl;
        exit(-1);
    }

    int int_size = sizeof(int);
    size_t vertex_array_bytes = ((size_t)vertices_count_) * 4;
    deg_outputfile.write(reinterpret_cast<const char *>(&int_size), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&vertices_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(&edges_count_), 4);
    deg_outputfile.write(reinterpret_cast<const char *>(degrees), vertex_array_bytes);

    deg_outputfile.close();
    deg_outputfile.clear();

    delete[] degrees;

    std::ofstream edge_outputfile(edge_path, std::ios::binary);

    if (edge_outputfile.is_open()) {
        std::cout << "Open edge file " << edge_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot edge file " << edge_path << " ." << std::endl;
        exit(-1);
    }

    size_t edge_array_bytes = ((size_t)edges_count_ * 2) * 4;
    edge_outputfile.write(reinterpret_cast<const char *>(neighbors_), edge_array_bytes);

    edge_outputfile.close();
    edge_outputfile.clear();

    if (edge_labels_ != nullptr && !edge_label_file_path.empty()) {
        std::ofstream edge_label_outputfile(edge_label_file_path, std::ios::binary);

        if (edge_label_outputfile.is_open()) {
            std::cout << "Open edge label file " << edge_label_file_path << " successfully." << std::endl;
        } else {
            std::cerr << "Cannot open edge label file " << edge_label_file_path << " ." << std::endl;
            exit(-1);
        }

        edge_label_outputfile.write(reinterpret_cast<const char *>(edge_labels_), edge_array_bytes);

        edge_label_outputfile.close();
        edge_label_outputfile.clear();
    }

    std::ofstream label_outputfile(label_path, std::ios::binary);

    if (label_outputfile.is_open()) {
        std::cout << "Open label file " << label_path << " successfully." << std::endl;
    }
    else {
        std::cerr << "Cannot label file " << label_path << " ." << std::endl;
        exit(-1);
    }

    size_t label_array_bytes = ((size_t)vertices_count_) * 4;
    label_outputfile.write(reinterpret_cast<const char *>(vertex_labels_), label_array_bytes);

    label_outputfile.close();
    label_outputfile.clear();
}

void Graph::buildEdgeIndex() {
    edge_index_ = new sparse_hash_map<uint64_t, std::vector<edge>*>();

    edge cur_edge;
    for (uint32_t u = 0; u < vertices_count_; ++u) {
        uint32_t u_l = getVertexLabel(u);

        uint32_t u_nbrs_cnt;
        const uint32_t* u_nbrs = getVertexNeighbors(u, u_nbrs_cnt);

        cur_edge.vertices_[0] = u;

        for (uint32_t i = 0; i < u_nbrs_cnt; ++i) {
            uint32_t v = u_nbrs[i];
            uint32_t v_l = getVertexLabel(v);

            uint64_t key = (uint64_t) u_l << 32 | v_l;
            cur_edge.vertices_[1] = v;

            if (!edge_index_->contains(key)) {
                (*edge_index_)[key] = new std::vector<edge>();
            }
            (*edge_index_)[key]->emplace_back(cur_edge);
        }
    }

}

void Graph::loadGraphFromFileCompressed(const std::string &file_path) {
    std::string degree_file_path = file_path + "/" + degree_file_name;
    std::string adj_file_path = file_path + "/" + adj_file_name;
    std::string vertex_label_file_path = file_path + "/" + vertex_label_file_name;
    std::string edge_label_file_path = file_path + "/" + edge_label_file_name;

    loadGraphFromFileCompressed(degree_file_path, adj_file_path, vertex_label_file_path, edge_label_file_path);
}

void Graph::storeCompressedGraph(const std::string &file_path) {
    std::string degree_file_path = file_path + "/" + degree_file_name;
    std::string adj_file_path = file_path + "/" + adj_file_name;
    std::string vertex_label_file_path = file_path + "/" + vertex_label_file_name;
    std::string edge_label_file_path = file_path + "/" + edge_label_file_name;

    storeCompressedGraph(degree_file_path, adj_file_path, vertex_label_file_path, edge_label_file_path);
}

void Graph::loadGraphFromFileWithoutMeta(const std::string &file_path) {
    typedef struct TempEdge {
        VertexID begin_;
        VertexID end_;
        VertexID edge_label_;
    } TempEdge;
    spp::sparse_hash_map<VertexID, LabelID> vertices;
    vertices.reserve(4 * 1024 * 1024);
    std::vector<TempEdge> edges;
    edges.reserve(128 * 1024 * 1024);

    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    VertexID max_vertex_id = 0;
    LabelID max_label_id = 0;
    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID label;
            infile >> id >> label;
            vertices[id] = label;

            if (id > max_vertex_id) {
                max_vertex_id = id;
            }

            if (labels_frequency_.find(label) == labels_frequency_.end()) {
                labels_frequency_[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency_[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            TempEdge temp_edge;

            if (is_edge_labeled) {
                infile >> temp_edge.begin_ >> temp_edge.end_ >> temp_edge.edge_label_;
            }
            else {
                infile >> temp_edge.begin_ >> temp_edge.end_;
            }

            if (temp_edge.begin_ == temp_edge.end_) {
                continue;
            }
            else if (temp_edge.begin_ > temp_edge.end_) {
                std::swap(temp_edge.begin_, temp_edge.end_);
            }

            edges.emplace_back(temp_edge);
        }
    }

    infile.close();

    labels_count_ =
            (ui) labels_frequency_.size() > (max_label_id + 2) ? (ui) labels_frequency_.size() : max_label_id + 2;

    max_label_frequency_ = 0;
    for (auto element : labels_frequency_) {
        if (element.second > max_label_frequency_) {
            max_label_frequency_ = element.second;
        }
    }

    {
        // Convert to CSR.
        vertices_count_ = max_vertex_id + 1;
        vertex_labels_ = new LabelID[vertices_count_];
        std::fill(vertex_labels_, vertex_labels_ + vertices_count_, max_label_id + 1);
        for (auto& kv : vertices) {
            vertex_labels_[kv.first] = kv.second;
        }

        std::vector<ui> temp_degree(vertices_count_);
        std::sort(edges.begin(), edges.end(), [](const TempEdge& l, const TempEdge& r) {
           if (l.begin_ == r.begin_)
               return l.end_ < r.end_;
           return l.begin_ < r.begin_;
        });

        ui prev_begin = 0;
        ui prev_end = 0;

        for (auto& edge : edges) {
            // Remove duplicate edge.
            if (edge.begin_ == prev_begin && edge.end_ == prev_end) {
                continue;
            }
            edges_count_ += 1;

            temp_degree[edge.begin_] += 1;
            temp_degree[edge.end_] += 1;

            prev_begin = edge.begin_;
            prev_end = edge.end_;
        }

        offsets_ = new ui[vertices_count_ + 1];
        offsets_[0] = 0;
        for (ui i = 0; i < vertices_count_; ++i) {
            if (temp_degree[i] > max_degree_) {
                max_degree_ = temp_degree[i];
            }
            offsets_[i + 1] = offsets_[i] + temp_degree[i];
        }

        neighbors_ = new ui[edges_count_ * 2];
        if (is_edge_labeled)
            edge_labels_ = new LabelID[edges_count_ * 2];

        prev_begin = 0;
        prev_end = 0;
        std::fill(temp_degree.begin(), temp_degree.end(), 0);

        for (auto& edge : edges) {
            // Remove duplicate edge.
            if (edge.begin_ == prev_begin && edge.end_ == prev_end) {
                continue;
            }
            prev_begin = edge.begin_;
            prev_end = edge.end_;

            ui begin_offset = offsets_[edge.begin_] + temp_degree[edge.begin_];
            ui end_offset = offsets_[edge.end_] + temp_degree[edge.end_];

            neighbors_[begin_offset] = edge.end_;
            neighbors_[end_offset] = edge.begin_;

            if (is_edge_labeled) {
                edge_labels_[begin_offset] = edge.edge_label_;
                edge_labels_[end_offset] = edge.edge_label_;
            }

            temp_degree[edge.begin_] += 1;
            temp_degree[edge.end_] += 1;
        }
    }

    BuildReverseIndex();
    buildEdgeIndex();

#if OPTIMIZED_LABELED_GRAPH == 1
        if (enable_label_offset_) {
            BuildNLF();
            BuildLabelOffset();
        }
#endif
}
