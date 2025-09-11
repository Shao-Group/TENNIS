"""
BSD 3-Clause License

Copyright (c) 2025, Xiaofei Carl Zang, Ke Chen, Mingfu Shao, and The Pennsylvania State University

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import networkx as nx
import numpy as np
from typing import List, Union


class boundComputer:
    def __init__(self, matrix: Union[List[List[int]], np.ndarray] = None):
        self.matrix = matrix
        self.graph = None
        self.num_connected_components = 0
        self.hamming_distances = None  # Store all pairwise hamming distances
        self.connected_components = None  # Store connected components
        self.hypergraph = None  # Hypergraph where nodes are CCs
        self.upper_bound_mst = None  # Upper bound using MST approach
        self.upper_bound_hub = None  # Upper bound using hub (OR) approach
        self.upper_bound = None  # Final upper bound (minimum of both)

        self.sanity_check()
        if matrix is not None:
            self._comp_hamming_distance()
            self._build_graph()
            self._count_connected_components()
            self._build_hypergraph()
            self._compute_upper_bound_mst()
            self._compute_upper_bound_hub()
            self._compute_final_upper_bound()

    def sanity_check(self):
        """All vectors must be binary and of the same length. Convert to numpy array."""
        assert self.matrix is not None
        # Convert to numpy array for efficient operations
        self.matrix = np.array(self.matrix)
        # Check that all elements are binary (0 or 1)
        assert np.all(np.isin(self.matrix, [0, 1])), "All elements must be 0 or 1"
        assert self.matrix.ndim == 2, "Input must be a 2D array"
        return 0
    
    def _comp_hamming_distance(self):
        """Compute hamming distances for all pairs and store them."""
        matrix_expanded1 = self.matrix[:, np.newaxis, :]
        matrix_expanded2 = self.matrix[np.newaxis, :, :]
        self.hamming_distances = np.sum(matrix_expanded1 != matrix_expanded2, axis=2)
        return 0

    def _build_graph(self) -> None:
        """Build graph where nodes are vectors and edges connect vectors differing by 1 bit."""
        self.graph = nx.Graph()
        
        n_vectors = len(self.matrix)
        
        # Add nodes (each vector becomes a node, using tuple for hashability)
        for i, vector in enumerate(self.matrix):
            self.graph.add_node(i, vector=tuple(vector))
        
        # Find pairs with exactly distance 1 and add edges
        upper_tri_mask = np.zeros((n_vectors, n_vectors), dtype=bool)
        upper_tri_i, upper_tri_j = np.triu_indices(n_vectors, k=1)
        upper_tri_mask[upper_tri_i, upper_tri_j] = True

        i_indices, j_indices = np.where((self.hamming_distances == 1) & upper_tri_mask)
        
        for i, j in zip(i_indices, j_indices):
            self.graph.add_edge(i, j)
    
    def _count_connected_components(self) -> None:
        """Count the number of connected components in the graph."""
        if self.graph is None:
            raise ValueError("Graph is not built")
        
        self.num_connected_components = nx.number_connected_components(self.graph)
        self.connected_components = [list(component) for component in nx.connected_components(self.graph)]
    
    def get_connected_components_count(self) -> int:
        """Return the number of connected components."""
        return self.num_connected_components
    
    def get_connected_components(self) -> List[List[int]]:
        """Return list of connected components (each component is a list of node indices)."""
        return self.connected_components
    
    def _build_hypergraph(self) -> None:
        """Build hypergraph where nodes are connected components and edges have min distance between CCs."""
        if self.connected_components is None or self.hamming_distances is None:
            raise ValueError("Connected components and hamming distances must be computed first")
        n_components = len(self.connected_components)
        self.hypergraph = nx.Graph()
        
        # Add nodes (each connected component becomes a node)
        for i, component in enumerate(self.connected_components):
            self.hypergraph.add_node(i, component=component, size=len(component))
        
        # Add edges between all pairs of components with minimum hamming distance as weight
        for i in range(n_components):
            for j in range(i + 1, n_components):
                min_distance = self._min_distance_between_components(
                    self.connected_components[i], 
                    self.connected_components[j]
                )
                self.hypergraph.add_edge(i, j, weight=min_distance)
    
    def _min_distance_between_components(self, component1: List[int], component2: List[int]) -> int:
        """Find minimum hamming distance between any two vectors from different components."""
        min_dist = float('inf')
        
        for node1 in component1:
            for node2 in component2:
                # Use stored hamming distances (symmetric matrix)
                dist = self.hamming_distances[node1, node2]
                min_dist = min(min_dist, dist)
        
        return int(min_dist)
    
    def _compute_upper_bound_mst(self) -> None:
        """Compute upper bound based on hypergraph MST structure."""
        if self.hypergraph is None:
            raise ValueError("Hypergraph must be built first")
        
        # If there's only one connected component, upper bound is 0
        if self.num_connected_components <= 1:
            self.upper_bound_mst = 0
            return
        
        # Find minimum spanning tree of hypergraph to get minimal connection cost
        mst = nx.minimum_spanning_tree(self.hypergraph, weight='weight')
        
        # MST size: exactly (num_connected_components - 1) edges
        mst_size = len(mst.edges())
        assert mst_size == self.num_connected_components - 1, f"MST should have {self.num_connected_components - 1} edges, got {mst_size}"
        
        # Upper bound: sum of (min_distance - 1) for each MST edge
        # Each edge with weight d requires (d-1) intermediate nodes to connect components
        self.upper_bound_mst = sum(max(0, data['weight'] - 1) for _, _, data in mst.edges(data=True))
    
    def _compute_upper_bound_hub(self) -> None:
        """Compute upper bound using hub approach (bitwise OR + minimum flips per CC)."""
        if self.connected_components is None:
            raise ValueError("Connected components must be computed first")
        
        # If there's only one connected component, upper bound is 0
        if self.num_connected_components <= 1:
            self.upper_bound_hub = 0
            return
        
        # Compute bitwise OR of all vectors (positions with at least one 1)
        hub_vector = np.any(self.matrix, axis=0).astype(int)
        
        total_flips = 0
        
        # For each connected component, find the minimum flips needed to reach hub_vector
        for component in self.connected_components:
            min_flips_in_component = float('inf')
            
            for node_idx in component:
                # Calculate hamming distance from this vector to hub_vector
                flips_needed = np.sum(self.matrix[node_idx] != hub_vector)
                min_flips_in_component = min(min_flips_in_component, flips_needed)
            
            total_flips += min_flips_in_component
        
        self.upper_bound_hub = total_flips
    
    def _compute_final_upper_bound(self) -> None:
        """Compute final upper bound as minimum of MST and hub approaches."""
        if self.upper_bound_mst is None or self.upper_bound_hub is None:
            raise ValueError("Both MST and hub upper bounds must be computed first")
        
        # Take the minimum of both approaches
        self.upper_bound = min(self.upper_bound_mst, self.upper_bound_hub)
    
    def get_upper_bound(self) -> int:
        """Return the computed upper bound (minimum of MST and hub approaches)."""
        return self.upper_bound
    
    def get_upper_bound_mst(self) -> int:
        """Return the MST-based upper bound."""
        return self.upper_bound_mst
    
    def get_upper_bound_hub(self) -> int:
        """Return the hub-based upper bound."""
        return self.upper_bound_hub
    
    def get_hub_vector(self) -> np.ndarray:
        """Return the hub vector (bitwise OR of all input vectors)."""
        return np.any(self.matrix, axis=0).astype(int)
    
    def get_hypergraph(self) -> nx.Graph:
        """Return the hypergraph where nodes are connected components."""
        return self.hypergraph
    
    def get_mst_info(self) -> dict:
        """Return MST information including size and edge details."""
        if self.hypergraph is None:
            raise ValueError("Hypergraph must be built first")
        
        if self.num_connected_components <= 1:
            return {"mst_size": 0, "mst_edges": [], "total_weight": 0, "upper_bound_contribution": 0}
        
        mst = nx.minimum_spanning_tree(self.hypergraph, weight='weight')
        mst_edges = [(i, j, data['weight']) for i, j, data in mst.edges(data=True)]
        total_weight = sum(data['weight'] for _, _, data in mst.edges(data=True))
        
        return {
            "mst_size": len(mst.edges()),
            "mst_edges": mst_edges,
            "total_weight": total_weight,
            "upper_bound_contribution": sum(max(0, weight - 1) for _, _, weight in mst_edges)
        }
    
    def set_matrix(self, matrix: Union[List[List[int]], np.ndarray]) -> None:
        """Set a new matrix and rebuild the graph."""
        self.matrix = matrix
        self.sanity_check()
        self._build_graph()
        self._count_connected_components()
        self._build_hypergraph()
        self._compute_upper_bound_mst()
        self._compute_upper_bound_hub()
        self._compute_final_upper_bound()


def test_bound_computer():
    """Test function to demonstrate boundComputer functionality."""
    # Example binary matrix
    test_matrix = [
        [1, 0, 0],  # Node 0
        [1, 1, 0],  # Node 1 - differs from 0 by 1 bit
        [0, 1, 0],  # Node 2 - differs from 1 by 1 bit
        [0, 0, 1],  # Node 3 - isolated (differs from all others by >1 bit)
        [1, 0, 1]   # Node 4 - differs from 0 by 1 bit
    ]
    
    bc = boundComputer(test_matrix)
    
    print("Matrix:")
    for i, row in enumerate(test_matrix):
        print(f"  Node {i}: {row}")
    
    print(f"\nNumber of connected components: {bc.get_connected_components_count()}")
    print(f"Connected components: {bc.get_connected_components()}")
    
    print("\nHypergraph edges (component pairs with minimum distances):")
    for i, j, data in bc.get_hypergraph().edges(data=True):
        print(f"  Component {i} ↔ Component {j}: min distance = {data['weight']}")
    
    print("\nMST Information:")
    mst_info = bc.get_mst_info()
    print(f"  MST size (number of edges): {mst_info['mst_size']}")
    print(f"  MST edges: {mst_info['mst_edges']}")
    print(f"  Total MST weight: {mst_info['total_weight']}")
    print(f"  Upper bound contribution: {mst_info['upper_bound_contribution']}")
    
    print(f"\nHub Vector (bitwise OR): {bc.get_hub_vector()}")
    
    print("\nUpper Bound Analysis:")
    print(f"  MST-based upper bound: {bc.get_upper_bound_mst()}")
    print(f"  Hub-based upper bound: {bc.get_upper_bound_hub()}")
    print(f"  Final upper bound (minimum): {bc.get_upper_bound()}")
    
    # Expected result: 2 connected components
    # Component 1: [0, 1, 2, 4] (nodes 0-1-2 connected, and 0-4 connected)  
    # Component 2: [3] (node 3 is isolated)
    # Hub vector: [1, 1, 1] (OR of all vectors)
    # MST approach: minimum distance between components - 1
    # Hub approach: sum of minimum flips from each CC to hub vector


if __name__ == "__main__":
    test_bound_computer()
        