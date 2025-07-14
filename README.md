# RICH

## Introduction

Arbitrage is a challenging data science problem characterized by rapidly fluctuating price discrepancies across token pairs, necessitating real-time solutions. To address this challenge, we model it as a \(k\)-hop negative cycle detection problem in token graphs and introduce **RICH**: **R**eal-time **I**dentification of negative **C**ycles for **H**igh-efficiency arbitrage.


## Repository Overview

This repository contains the implementation of **RICH**, including the source code and sample datasets. The code consists of:

- A preprocessing script for graph conversion.
- A compiled executable for negative cycle detection.
- Scripts for running experiments on sample datasets.

## Compilation and Execution

### 1. Convert Graph to CSR Format

Before running **RICH**, the input graph must be converted into a compressed sparse row (CSR) format. Run the following command:

```zsh
python convert_to_csr.py original_graph_path csr_graph_path
```

### 2. Compile the Code

Compile the source code by running:

```zsh
mkdir build
cd build
cmake ..
make
```

### 3. Run the Executable

After compilation, run the binary for negative cycle detection:

```zsh
./DPWithFilterFull
```

If the at-most k variant is preferred, then use:

```zsh
./DPWithFilterFullAtMostk
```

Then, specify the dataset, k, number of coloring process. 

## Data

This repository includes sample datasets to test **RICH**. The datasets contain transaction-based graphs where nodes represent tokens and edges represent token exchange interactions.

## Usage Example

To execute **RICH** on a dataset, follow these steps:

1. Convert the input graph into CSR format.
2. Compile the source code.
3. Run the compiled executable on the processed graph.


## Contact

For questions or collaboration inquiries, feel free to contact us.


## Disclaimer

RICH is developed exclusively for academic research. Its deployment in real-world trading scenarios must account for regulatory and ethical considerations. Any practical use is at the user's own risk and should carefully consider multiple market factors such as slippage and execution delays.

---

