#include "directed_graph.h"
#include <limits>
#include <iostream>
#include <random>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <string>
#include <fstream> // Include fstream for file output
#include <sstream>
#include <queue>

// Custom hash function for std::vector<int>
struct VectorHash {
    std::size_t operator()(const std::vector<int>& vec) const {
        std::size_t seed = vec.size();
        for (auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// KCycleColorCoding class for detecting negative-weight k-cycles
class KCycleColorCoding {
public:
    KCycleColorCoding(const DirectedGraph* graph, int k, int num_trials, unsigned int seed)
        : graph_(graph), k_(k), num_trials_(num_trials), gen_(seed) {}

    std::pair<double, std::vector<uint32_t>> find_most_negative_k_cycle(
        unsigned int seed, 
        std::ofstream& results_file,
        const std::vector<unsigned int>& trial_seeds);

private:
    const DirectedGraph* graph_;
    DirectedGraph new_graph;
    int k_;
    int num_trials_;
    std::mt19937 gen_;

    std::unordered_map<uint32_t, int> assign_colors(const std::vector<uint32_t>& vertices);

    void bfs_dp_search_k_cycle(
        uint32_t start,
        // double current_weight,
        // int depth,
        // std::vector<uint32_t>& path,
        // int color_set,
        const std::unordered_map<uint32_t,int>& colors,
        double& local_best_weight,
        std::vector<uint32_t>& local_best_cycle,
        std::vector<bool>& fully_explored,
        std::unordered_map<uint32_t, std::unordered_map<int, double>>& dp_table
    );
};

std::unordered_map<uint32_t, int> KCycleColorCoding::assign_colors(const std::vector<uint32_t>& vertices) {
    std::unordered_map<uint32_t, int> colors;
    std::uniform_int_distribution<int> dist(0, k_ - 1);
    // for (auto v : vertices) {
    for(uint32_t v = 0; v < vertices.size(); v++) {
        colors[v] = dist(gen_);
        // colors[v] = v % k_;
    }
    return colors;
}


void KCycleColorCoding::bfs_dp_search_k_cycle(
    uint32_t start,
    const std::unordered_map<uint32_t,int>& colors,
    double& local_best_weight,
    std::vector<uint32_t>& local_best_cycle,
    std::vector<bool>& fully_explored,
    std::unordered_map<uint32_t, std::unordered_map<int, double>>& dp_table
) {
    std::queue<std::tuple<uint32_t, int, double, std::vector<uint32_t>>> q;

    int start_color = 1 << colors.at(start);
    dp_table[start][start_color] = 0.0;
    
    std::vector<uint32_t> init_path = {start};
    q.push({start, start_color, 0.0, init_path});

    while (!q.empty()) {
        auto [current, color_set, current_weight, path] = q.front();
        q.pop();

        if (path.size() == k_) {
            auto [neighbors, weights, deg] = new_graph.get_out_neighbors(current);
            for (int i = 0; i < deg; i++) {
                if (neighbors[i] == start) { 
                    double cycle_weight = current_weight + weights[i];

                    if (cycle_weight < local_best_weight) {
                        local_best_weight = cycle_weight;
                        local_best_cycle = path;
                        local_best_cycle.push_back(start);
                    }
                }
            }
            continue; 
        }

        auto [neighbors, weights, deg] = new_graph.get_out_neighbors(current);
        for (int i = 0; i < deg; i++) {
            uint32_t nxt = neighbors[i];
            double w = weights[i];

            if (fully_explored[nxt]) {
                continue;  
            }

            int nxt_color = 1 << colors.at(nxt);
            if (color_set & nxt_color) {
                continue; 
            }

            int new_color_set = color_set | nxt_color;
            double new_weight = current_weight + w;

            if (!dp_table[nxt].count(new_color_set) || dp_table[nxt][new_color_set] > new_weight) {
                dp_table[nxt][new_color_set] = new_weight;

                std::vector<uint32_t> new_path = path;
                new_path.push_back(nxt);

                q.push({nxt, new_color_set, new_weight, new_path});
            }
        }
    }

    fully_explored[start] = true;
}


std::pair<double, std::vector<uint32_t>> KCycleColorCoding::find_most_negative_k_cycle(
    unsigned int seed, 
    std::ofstream& results_file,
    const std::vector<unsigned int>& trial_seeds) {
    // Prepare a list of all vertices.
    std::vector<uint32_t> vertices(graph_->num_vertices());
    for (uint32_t v = 0; v < graph_->num_vertices(); v++) {
        vertices[v] = v;
    }

    double global_most_negative_weight = std::numeric_limits<double>::infinity();
    std::vector<uint32_t> global_most_negative_cycle;

    // Run multiple trials with predetermined seeds
    for (int trial = 0; trial < num_trials_; ++trial) {
        unsigned int trial_seed = trial_seeds[trial];
        std::cout << "trial_seed: " << trial_seed << std::endl;
        gen_.seed(trial_seed);
        auto start_time = std::chrono::high_resolution_clock::now();

        // Assign colors randomly to each vertex in [0..k_-1].
        std::unordered_map<uint32_t, int> colors = assign_colors(vertices);
        // prune same color neighbor
        // Prune nodes based on color and neighbors.
        bool has_pruned;
        std::vector<bool> to_prune(vertices.size(), false);
        do {
            has_pruned = false;
            for (auto it = vertices.begin(); it != vertices.end(); ++it) {
                uint32_t a = *it;
                if (to_prune[a]) continue;

                auto [neighbors_ptr, weights, deg] = graph_->get_out_neighbors(a);
                std::vector<uint32_t> valid_neighbors;
                for (uint32_t i = 0; i < deg; i++) {
                    if (!to_prune[neighbors_ptr[i]]) {
                        valid_neighbors.push_back(neighbors_ptr[i]);
                    }
                }
                if (valid_neighbors.size() == 1 && colors[a] == colors[valid_neighbors[0]]) {
                    to_prune[a] = true;
                    has_pruned = true;
                }
            }
        } while (has_pruned);
        std::cout << "pruned nodes: " << std::count(to_prune.begin(), to_prune.end(), true) << std::endl;
        // Create a new graph in memory by filtering out pruned nodes.
        //;
        new_graph = DirectedGraph();
        new_graph.num_vertices_ = std::count(to_prune.begin(), to_prune.end(), false);
        new_graph.num_edges_ = 0;
        std::unordered_map<uint32_t, uint32_t> old_to_new;
        std::unordered_map<uint32_t, uint32_t> new_to_old;
        uint32_t new_idx = 0;

        for (uint32_t i = 0; i < graph_->num_vertices(); ++i) {
            if (!to_prune[i]) {
                old_to_new[i] = new_idx;
                new_to_old[new_idx] = i;
                new_idx++;
            }
        }

        std::vector<uint32_t> new_out_degree(new_graph.num_vertices_, 0);
        std::vector<uint32_t> new_in_degree(new_graph.num_vertices_, 0);

        for (uint32_t u = 0; u < graph_->num_vertices(); ++u) {
            if (to_prune[u]) continue;

            auto [neighbors_ptr, weights, deg] = graph_->get_out_neighbors(u);
            for (uint32_t i = 0; i < deg; ++i) {
                uint32_t v = neighbors_ptr[i];
                if (!to_prune[v]) {
                    ++new_out_degree[old_to_new[u]];
                    ++new_in_degree[old_to_new[v]];
                    ++new_graph.num_edges_;
                }
            }
        }

        new_graph.out_degree_ = new uint32_t[new_graph.num_vertices_];
        std::copy(new_out_degree.begin(), new_out_degree.end(), new_graph.out_degree_);

        new_graph.in_degree_ = new uint32_t[new_graph.num_vertices_];
        std::copy(new_in_degree.begin(), new_in_degree.end(), new_graph.in_degree_);

        new_graph.out_offset_ = new uint32_t[new_graph.num_vertices_ + 1]();
        new_graph.in_offset_ = new uint32_t[new_graph.num_vertices_ + 1]();

        for (uint32_t i = 1; i <= new_graph.num_vertices_; ++i) {
            new_graph.out_offset_[i] = new_graph.out_offset_[i - 1] + new_out_degree[i - 1];
            new_graph.in_offset_[i] = new_graph.in_offset_[i - 1] + new_in_degree[i - 1];
        }

        uint32_t total_edges = new_graph.out_offset_[new_graph.num_vertices_];
        new_graph.out_adj_ = new uint32_t[total_edges];
        new_graph.out_weights_ = new double[total_edges];

        total_edges = new_graph.in_offset_[new_graph.num_vertices_];
        new_graph.in_adj_ = new uint32_t[total_edges];
        new_graph.in_weights_ = new double[total_edges];

        std::vector<uint32_t> out_counter(new_graph.out_offset_, new_graph.out_offset_ + new_graph.num_vertices_);
        std::vector<uint32_t> in_counter(new_graph.in_offset_, new_graph.in_offset_ + new_graph.num_vertices_);

        for (uint32_t u = 0; u < graph_->num_vertices(); ++u) {
            if (to_prune[u]) continue;

            auto [neighbors_ptr, weights, deg] = graph_->get_out_neighbors(u);
            for (uint32_t i = 0; i < deg; ++i) {
                uint32_t v = neighbors_ptr[i];
                if (!to_prune[v]) {
                    new_graph.out_adj_[out_counter[old_to_new[u]]] = old_to_new[v];
                    new_graph.out_weights_[out_counter[old_to_new[u]]++] = weights[i];

                    new_graph.in_adj_[in_counter[old_to_new[v]]] = old_to_new[u];
                    new_graph.in_weights_[in_counter[old_to_new[v]]++] = weights[i];
                }
            }
        }

        // Transfer colors to the new graph
        std::unordered_map<uint32_t, int> new_colors;
        for (const auto& [old_idx, new_idx] : old_to_new) {
            new_colors[new_idx] = colors[old_idx];
        }

        std::vector<uint32_t> new_vertices(new_graph.num_vertices());
        for (uint32_t i = 0; i < new_graph.num_vertices(); ++i) {
            new_vertices[i] = i;
        }

        std::cout << "new_graph_size: " << new_graph.num_vertices_ << std::endl;

        auto color_assign_time = std::chrono::high_resolution_clock::now();

        
        double trial_best_weight = std::numeric_limits<double>::infinity();
        std::vector<uint32_t> trial_best_cycle;

        std::vector<bool> fully_explored(new_graph.num_vertices(), false);

        for (auto start_node : new_vertices) {
            std::unordered_map<uint32_t, std::unordered_map<int, double>> dp_table;
            if (fully_explored[start_node]) {
                continue;
            }

            int color_set = 1 << new_colors[start_node];

            std::vector<uint32_t> path;
            path.push_back(start_node);

            // DP
            bfs_dp_search_k_cycle(
                start_node, new_colors,
                trial_best_weight, trial_best_cycle,
                fully_explored, dp_table
            );
            fully_explored[start_node] = true;
        }
        if (trial_best_weight < global_most_negative_weight) {
            global_most_negative_weight = trial_best_weight;
            global_most_negative_cycle.clear();
            for (uint32_t v : trial_best_cycle) {
                global_most_negative_cycle.push_back(new_to_old[v]);
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                                  end_time - start_time).count();
        auto color_assign_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                                         color_assign_time - start_time).count();

        results_file << "Trial " << (trial + 1) << ":\n";
        if (!global_most_negative_cycle.empty()
            && global_most_negative_weight < std::numeric_limits<double>::infinity()) 
        {
            results_file << "  Found a candidate cycle of weight: " 
                      << global_most_negative_weight << "\n  Cycle: ";
            for (size_t i = 0; i < global_most_negative_cycle.size(); ++i) {
                results_file << global_most_negative_cycle[i];
                if (i < global_most_negative_cycle.size() - 1)
                    results_file << " -> ";
            }
            results_file << "\n";
        } else {
            results_file << "  No valid " << k_ << "-cycle found.\n";
        }
        results_file << "  Color assignment time: " << color_assign_duration << " ms\n";
        results_file << "  Total trial time: " << total_duration << " ms\n";
        results_file << "  Actual trial time: " << total_duration - color_assign_duration << " ms\n";
    }

    return {global_most_negative_weight, global_most_negative_cycle};
}


int main() {
    std::cout << "Enter the graph directory: ";
    std::string graph_dir;
    std::cin >> graph_dir;

    DirectedGraph graph;
    graph.load_csr(graph_dir);
    graph.print_metadata();

    int k;
    std::cout << "\nEnter k (cycle length to detect): ";
    std::cin >> k;

    int num_trials;
    std::cout << "Enter number of trials: ";
    std::cin >> num_trials;

    // Generate random seeds
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> dist(0, std::numeric_limits<unsigned int>::max());
    
    // Generate one seed for the experiment
    unsigned int experiment_seed = dist(gen);
    
    // Generate seeds for trials
    std::vector<unsigned int> trial_seeds;
    for (int i = 0; i < num_trials; i++) {
        trial_seeds.push_back(dist(gen));
    }

    std::string results_filename = "../results/dp_full/" + std::to_string(k) + "_" + std::to_string(num_trials) + ".txt";
    std::ofstream results_file(results_filename);

    std::cout << "Running experiment with seed: " << experiment_seed << std::endl;

    KCycleColorCoding detector(&graph, k, num_trials, experiment_seed);
    auto result = detector.find_most_negative_k_cycle(experiment_seed, results_file, trial_seeds);

    double weight = result.first;
    std::vector<uint32_t> cycle = result.second;

    if (!cycle.empty()) {
        results_file << "Most negative cycle weight: " << weight << "\nCycle: ";
        for (size_t i = 0; i < cycle.size(); ++i) {
            results_file << cycle[i];
            if (i < cycle.size() - 1)
                results_file << " -> ";
        }
        results_file << "\n";
    } else {
        results_file << "No negative cycle found.\n";
    }

    results_file.close();
    return 0;
}