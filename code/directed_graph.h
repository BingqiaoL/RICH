#ifndef DIRECTED_CSR_H
#define DIRECTED_CSR_H

#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <cstdint>

class DirectedGraph {
public:
    // Meta information
    uint32_t num_vertices_;
    uint32_t num_edges_;
    uint32_t max_num_in_neighbors_;
    uint32_t max_num_out_neighbors_;

    // CSR representation
    uint32_t* out_offset_;
    uint32_t* out_adj_;
    double* out_weights_;
    uint32_t* out_degree_;

    uint32_t* in_offset_;
    uint32_t* in_adj_;
    double* in_weights_;
    uint32_t* in_degree_;

    // Constructor
    explicit DirectedGraph()
        : num_vertices_(0), num_edges_(0), max_num_in_neighbors_(0), max_num_out_neighbors_(0),
          out_offset_(nullptr), out_adj_(nullptr), out_weights_(nullptr), out_degree_(nullptr),
          in_offset_(nullptr), in_adj_(nullptr), in_weights_(nullptr), in_degree_(nullptr) {}

    // Destructor
    ~DirectedGraph() {
        free(out_offset_);
        free(in_offset_);
        free(out_adj_);
        free(out_weights_);
        free(out_degree_);
        free(in_adj_);
        free(in_weights_);
        free(in_degree_);
    }

    // Accessors
    inline uint32_t num_vertices() const {
        return num_vertices_;
    }

    inline uint32_t num_edges() const {
        return num_edges_;
    }

    /**
     * Get out-neighbors of vertex u along with their weights.
     * Returns a tuple containing a pointer to the adjacency list, a pointer to the weights, and the degree.
     */
    inline std::tuple<const uint32_t*, const double*, uint32_t> get_out_neighbors(uint32_t u) const {
        return std::make_tuple(out_adj_ + out_offset_[u], out_weights_ + out_offset_[u], out_degree_[u]);
    }

    /**
     * Get in-neighbors of vertex u along with their weights.
     * Returns a tuple containing a pointer to the adjacency list, a pointer to the weights, and the degree.
     */
    inline std::tuple<const uint32_t*, const double*, uint32_t> get_in_neighbors(uint32_t u) const {
        return std::make_tuple(in_adj_ + in_offset_[u], in_weights_ + in_offset_[u], in_degree_[u]);
    }

    inline uint32_t num_out_neighbors(uint32_t u) const {
        return out_degree_[u];
    }

    inline uint32_t num_in_neighbors(uint32_t u) const {
        return in_degree_[u];
    }

    inline uint32_t num_neighbors(uint32_t u) const {
        return in_degree_[u] + out_degree_[u];
    }

    // Graph I/O operations
    void load_csr(const std::string& graph_dir);
    void store_csr(const std::string& graph_dir);
    void load_edge_list(const std::string& graph_dir, char skip = '#');
    void store_edge_list(const std::string& graph_dir);
    void print_metadata();

private:
    // Graph I/O helper functions
    void load_csr_degree_file(const std::string& degree_file_path, uint32_t*& degree);
    void load_csr_adj_file(const std::string& adj_file_path, const uint32_t* degree, uint32_t*& offset,
                           uint32_t*& adj, double*& weights, const std::string& weight_file_path);
};

#endif // DIRECTED_CSR_H