#include "directed_graph.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cassert>

void DirectedGraph::load_csr(const std::string& graph_dir) {
    // Paths to CSR files
    std::string out_degree_file = graph_dir + "/b_out_degree.bin";
    std::string out_adj_file = graph_dir + "/b_out_adj.bin";
    std::string out_weight_file = graph_dir + "/b_out_weight.bin";

    std::string in_degree_file = graph_dir + "/b_in_degree.bin";
    std::string in_adj_file = graph_dir + "/b_in_adj.bin";
    std::string in_weight_file = graph_dir + "/b_in_weight.bin";

    // Load out-degree
    load_csr_degree_file(out_degree_file, out_degree_);

    // Determine number of vertices from out_degree_
    // Assuming that in_degree_ will have the same number of vertices
    // since the node IDs are remapped from 0 to max_num_nodes -1
    std::ifstream infile(out_degree_file, std::ios::binary | std::ios::ate);
    std::streamsize out_degree_size = infile.tellg();
    infile.close();
    num_vertices_ = out_degree_size / sizeof(uint32_t);

    // Update max_num_out_neighbors_
    max_num_out_neighbors_ = 0;
    for (uint32_t i = 0; i < num_vertices_; ++i) {
        if (out_degree_[i] > max_num_out_neighbors_) {
            max_num_out_neighbors_ = out_degree_[i];
        }
    }

    // Load out-adj and out-weights
    load_csr_adj_file(out_adj_file, out_degree_, out_offset_, out_adj_, out_weights_, out_weight_file);

    // Load in-degree
    load_csr_degree_file(in_degree_file, in_degree_);

    // Load in-adj and in-weights
    load_csr_adj_file(in_adj_file, in_degree_, in_offset_, in_adj_, in_weights_, in_weight_file);

    // Calculate number of edges and max_num_in_neighbors_
    num_edges_ = 0;
    max_num_in_neighbors_ = 0;
    for (uint32_t i = 0; i < num_vertices_; ++i) {
        num_edges_ += out_degree_[i];
        if (in_degree_[i] > max_num_in_neighbors_) {
            max_num_in_neighbors_ = in_degree_[i];
        }
    }
}

void DirectedGraph::load_csr_degree_file(const std::string& degree_file_path, uint32_t*& degree) {
    std::ifstream infile(degree_file_path, std::ios::binary | std::ios::ate);
    std::streamsize size = infile.tellg();
    infile.seekg(0, std::ios::beg);
    uint32_t num_vertices = size / sizeof(uint32_t);
    degree = (uint32_t*)malloc(num_vertices * sizeof(uint32_t));
    infile.read(reinterpret_cast<char*>(degree), size);
    infile.close();
}

void DirectedGraph::load_csr_adj_file(const std::string& adj_file_path, const uint32_t* degree, uint32_t*& offset,
                                     uint32_t*& adj, double*& weights, const std::string& weight_file_path) {
    // Allocate and compute offset
    offset = (uint32_t*)malloc((num_vertices_ + 1) * sizeof(uint32_t));
    offset[0] = 0;
    for (uint32_t i = 0; i < num_vertices_; ++i) {
        offset[i + 1] = offset[i] + degree[i];
    }

    uint32_t total_edges = offset[num_vertices_];
    adj = (uint32_t*)malloc(total_edges * sizeof(uint32_t));
    weights = (double*)malloc(total_edges * sizeof(double));

    // Load adjacency
    std::ifstream adj_infile(adj_file_path, std::ios::binary);
    adj_infile.read(reinterpret_cast<char*>(adj), total_edges * sizeof(uint32_t));
    adj_infile.close();

    // Load weights
    std::ifstream weight_infile(weight_file_path, std::ios::binary);
    weight_infile.read(reinterpret_cast<char*>(weights), total_edges * sizeof(double));
    weight_infile.close();
}

void DirectedGraph::print_metadata() {
    std::cout << "Graph Metadata:\n";
    std::cout << "Number of vertices: " << num_vertices_ << "\n";
    std::cout << "Number of edges: " << num_edges_ << "\n";
    std::cout << "Maximum out-degree: " << max_num_out_neighbors_ << "\n";
    std::cout << "Maximum in-degree: " << max_num_in_neighbors_ << "\n";
}

void DirectedGraph::store_csr(const std::string& graph_dir) {
    // Implementation can be added if needed
}

void DirectedGraph::load_edge_list(const std::string& graph_dir, char skip) {
    // Implementation can be added if needed
}

void DirectedGraph::store_edge_list(const std::string& graph_dir) {
    // Implementation can be added if needed
}