import argparse
import os
import numpy as np
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert edge list to CSR binary files with node remapping.')
    parser.add_argument('input_file', help='Path to the input edge list file.')
    parser.add_argument('output_folder', help='Path to the output folder.')
    return parser.parse_args()

def read_edges(input_file):
    edges = []
    original_nodes = set()
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Split by comma or whitespace
            if ',' in line:
                parts = line.split(',')
            else:
                parts = line.split()
            if len(parts) < 3:
                continue  # Skip invalid lines
            src, dst, weight = parts[:3]
            src = int(src)
            dst = int(dst)
            weight = float(weight)
            edges.append((src, dst, weight))
            original_nodes.add(src)
            original_nodes.add(dst)
    return edges, original_nodes

def create_node_mapping(original_nodes):
    sorted_nodes = sorted(original_nodes)
    node_map = {original_id: new_id for new_id, original_id in enumerate(sorted_nodes)}
    return node_map

def write_node_map(node_map, output_folder):
    map_file = os.path.join(output_folder, 'node_map.txt')
    with open(map_file, 'w') as f:
        for original_id, new_id in node_map.items():
            f.write(f"{original_id} {new_id}\n")

def remap_edges(edges, node_map):
    remapped = []
    for src, dst, weight in edges:
        new_src = node_map[src]
        new_dst = node_map[dst]
        remapped.append((new_src, new_dst, weight))
    return remapped

def build_csr(remapped_edges, num_vertices):
    # Sort edges by source
    remapped_edges.sort(key=lambda x: x[0])
    
    # Initialize out-degree and in-degree
    out_degree = np.zeros(num_vertices, dtype=np.uint32)
    in_degree = np.zeros(num_vertices, dtype=np.uint32)
    for src, dst, _ in remapped_edges:
        out_degree[src] += 1
        in_degree[dst] += 1

    # Build out_offset
    out_offset = np.zeros(num_vertices + 1, dtype=np.uint32)
    out_offset[1:] = np.cumsum(out_degree)

    # Build in_offset
    in_offset = np.zeros(num_vertices + 1, dtype=np.uint32)
    in_offset[1:] = np.cumsum(in_degree)

    # Initialize adjacency and weight arrays
    out_adj = np.zeros(len(remapped_edges), dtype=np.uint32)
    out_weights = np.zeros(len(remapped_edges), dtype=np.float64)

    in_adj = np.zeros(len(remapped_edges), dtype=np.uint32)
    in_weights = np.zeros(len(remapped_edges), dtype=np.float64)

    # Temporary counters to keep track of insertion points
    out_counter = out_offset[:-1].copy()
    in_counter = in_offset[:-1].copy()

    for src, dst, weight in remapped_edges:
        out_idx = out_counter[src]
        out_adj[out_idx] = dst
        out_weights[out_idx] = weight
        out_counter[src] += 1

        in_idx = in_counter[dst]
        in_adj[in_idx] = src
        in_weights[in_idx] = weight
        in_counter[dst] += 1

    return out_degree, out_offset, out_adj, out_weights, in_degree, in_offset, in_adj, in_weights

def write_binary_files(out_degree, out_adj, out_weights, in_degree, in_adj, in_weights, output_folder):
    out_degree.tofile(os.path.join(output_folder, 'b_out_degree.bin'))
    out_adj.tofile(os.path.join(output_folder, 'b_out_adj.bin'))
    out_weights.tofile(os.path.join(output_folder, 'b_out_weight.bin'))

    in_degree.tofile(os.path.join(output_folder, 'b_in_degree.bin'))
    in_adj.tofile(os.path.join(output_folder, 'b_in_adj.bin'))
    in_weights.tofile(os.path.join(output_folder, 'b_in_weight.bin'))

def main():
    args = parse_arguments()
    input_file = args.input_file
    output_folder = args.output_folder

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Step 1: Read edges and collect unique nodes
    edges, original_nodes = read_edges(input_file)

    # Step 2: Create node mapping
    node_map = create_node_mapping(original_nodes)

    # Step 3: Write node mapping to file
    write_node_map(node_map, output_folder)

    # Step 4: Remap edges with new node IDs
    remapped_edges = remap_edges(edges, node_map)

    # Step 5: Build CSR structures
    num_vertices = len(node_map)
    out_degree, out_offset, out_adj, out_weights, in_degree, in_offset, in_adj, in_weights = build_csr(remapped_edges, num_vertices)

    # Step 6: Write CSR binary files
    write_binary_files(out_degree, out_adj, out_weights, in_degree, in_adj, in_weights, output_folder)

    print("CSR files and node_map.txt have been successfully generated.")

if __name__ == "__main__":
    main()