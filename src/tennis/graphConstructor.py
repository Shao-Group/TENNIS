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

import numpy as np
import networkx as nx
from typing import List, Dict, Tuple, Union, Optional

try:
    from .util import *
except ImportError:
    # For standalone execution
    from util import *


class GraphConstructor:
    def __init__(self, original_nodes: List[List[int]], constructed_nodes: List[List[int]] = None):
        """
        Initialize graph constructor with original and newly constructed nodes.
        
        Args:
            original_nodes: List of original binary vectors
            constructed_nodes: List of newly constructed binary vectors (optional)
        """
        self.original_nodes = np.array(original_nodes)
        self.constructed_nodes = np.array(constructed_nodes) if constructed_nodes else np.empty((0, len(original_nodes[0])))
        
        # Combine all nodes
        if len(self.constructed_nodes) > 0:
            self.all_nodes = np.vstack([self.original_nodes, self.constructed_nodes])
        else:
            self.all_nodes = self.original_nodes.copy()
        
        self.num_original = len(self.original_nodes)
        self.num_constructed = len(self.constructed_nodes)
        self.num_total = len(self.all_nodes)
        
        # Build basic connectivity graph (hamming distance = 1)
        self.connectivity_graph = None
        self._build_connectivity_graph()
        
        # Weighted graph for evolutionary tree
        self.weighted_graph = None
        self.max_weight_tree = None
        
    def _build_connectivity_graph(self) -> None:
        """Build graph where nodes are connected if hamming distance = 1."""
        self.connectivity_graph = nx.Graph()
        
        # Add all nodes
        for i in range(self.num_total):
            node_type = "original" if i < self.num_original else "constructed"
            self.connectivity_graph.add_node(i, 
                                           vector=tuple(self.all_nodes[i]),
                                           type=node_type,
                                           index_in_type=i if node_type == "original" else i - self.num_original)
        
        # Add edges for nodes with hamming distance = 1
        for i in range(self.num_total):
            for j in range(i + 1, self.num_total):
                hamming_dist = np.sum(self.all_nodes[i] != self.all_nodes[j])
                if hamming_dist == 1:
                    self.connectivity_graph.add_edge(i, j, hamming_distance=1)
    
    def set_edge_weights(self, weight_function: callable) -> None:
        """
        Set edge weights using an external weight function.
        
        Args:
            weight_function: Function that takes (vector1, vector2, node_info1, node_info2) 
                           and returns a weight value
        """
        self.weighted_graph = self.connectivity_graph.copy()
        
        # Apply weights to all edges
        for i, j in self.weighted_graph.edges():
            vec1 = self.all_nodes[i]
            vec2 = self.all_nodes[j]
            info1 = self.weighted_graph.nodes[i]
            info2 = self.weighted_graph.nodes[j]
            
            weight = weight_function(vec1, vec2, info1, info2)
            self.weighted_graph[i][j]['weight'] = weight
    
    def set_edge_weights_from_dict(self, weight_dict: Dict[Tuple[int, int], float]) -> None:
        """
        Set edge weights from a dictionary mapping (node1, node2) -> weight.
        
        Args:
            weight_dict: Dictionary with (i, j) tuples as keys and weights as values
        """
        self.weighted_graph = self.connectivity_graph.copy()
        
        # Apply weights from dictionary
        for i, j in self.weighted_graph.edges():
            # Try both orientations since edges are undirected
            weight = weight_dict.get((i, j)) or weight_dict.get((j, i))
            if weight is not None:
                self.weighted_graph[i][j]['weight'] = weight
            else:
                # Default weight if not specified
                self.weighted_graph[i][j]['weight'] = 1.0
    
    def compute_maximum_weight_spanning_tree(self) -> nx.Graph:
        """
        Compute maximum weight spanning tree from the weighted graph.
        
        Returns:
            NetworkX graph representing the maximum weight spanning tree
        """
        if self.weighted_graph is None:
            raise ValueError("Edge weights must be set before computing spanning tree")
        
        # NetworkX minimum_spanning_tree finds minimum weight, so we negate weights
        negated_graph = self.weighted_graph.copy()
        for i, j in negated_graph.edges():
            negated_graph[i][j]['weight'] = -negated_graph[i][j]['weight']
        
        # Compute minimum spanning tree on negated weights = maximum on original weights
        mst = nx.minimum_spanning_tree(negated_graph, weight='weight')
        
        # Restore original weights
        for i, j in mst.edges():
            mst[i][j]['weight'] = -mst[i][j]['weight']
        
        self.max_weight_tree = mst
        return self.max_weight_tree
    
    def get_connectivity_graph(self) -> nx.Graph:
        """Return the basic connectivity graph (hamming distance = 1)."""
        return self.connectivity_graph
    
    def get_weighted_graph(self) -> Optional[nx.Graph]:
        """Return the weighted graph."""
        return self.weighted_graph
    
    def get_max_weight_tree(self) -> Optional[nx.Graph]:
        """Return the maximum weight spanning tree."""
        return self.max_weight_tree
    
    def get_tree_info(self) -> Dict:
        """Return information about the maximum weight spanning tree."""
        if self.max_weight_tree is None:
            return {"error": "Maximum weight tree not computed"}
        
        total_weight = sum(data['weight'] for _, _, data in self.max_weight_tree.edges(data=True))
        num_edges = len(self.max_weight_tree.edges())
        
        # Analyze tree structure
        original_nodes_in_tree = [n for n in self.max_weight_tree.nodes() if n < self.num_original]
        constructed_nodes_in_tree = [n for n in self.max_weight_tree.nodes() if n >= self.num_original]
        
        return {
            "total_weight": total_weight,
            "num_edges": num_edges,
            "num_nodes": len(self.max_weight_tree.nodes()),
            "original_nodes_in_tree": len(original_nodes_in_tree),
            "constructed_nodes_in_tree": len(constructed_nodes_in_tree),
            "is_connected": nx.is_connected(self.max_weight_tree),
            "tree_edges": [(i, j, data['weight']) for i, j, data in self.max_weight_tree.edges(data=True)]
        }
    
    def get_evolutionary_path(self, start_node: int, end_node: int) -> List[int]:
        """
        Find evolutionary path between two nodes in the maximum weight tree.
        
        Args:
            start_node: Starting node index
            end_node: Ending node index
            
        Returns:
            List of node indices representing the path
        """
        if self.max_weight_tree is None:
            raise ValueError("Maximum weight tree not computed")
        
        try:
            path = nx.shortest_path(self.max_weight_tree, start_node, end_node)
            return path
        except nx.NetworkXNoPath:
            return []
    
    def get_node_info(self, node_idx: int) -> Dict:
        """Return information about a specific node."""
        if node_idx >= self.num_total:
            raise ValueError(f"Node index {node_idx} out of range")
        
        node_type = "original" if node_idx < self.num_original else "constructed"
        type_index = node_idx if node_type == "original" else node_idx - self.num_original
        
        return {
            "index": node_idx,
            "type": node_type,
            "index_in_type": type_index,
            "vector": self.all_nodes[node_idx].tolist(),
            "neighbors_in_tree": list(self.max_weight_tree.neighbors(node_idx)) if self.max_weight_tree else []
        }


def test_graph_constructor():
    """Test function for the graph constructor."""
    # Original nodes (from connected components)
    original = [
        [1, 0, 0, 0],  # Node 0
        [1, 1, 0, 0],  # Node 1
        [0, 0, 1, 1],  # Node 2
        [0, 1, 1, 1],  # Node 3
    ]
    
    # Constructed intermediate nodes
    constructed = [
        [1, 0, 1, 0],  # Intermediate node 4
        [0, 1, 1, 0],  # Intermediate node 5
    ]
    
    gc = GraphConstructor(original, constructed)
    
    print("Graph Constructor Test")
    print(f"Original nodes: {len(gc.original_nodes)}")
    print(f"Constructed nodes: {len(gc.constructed_nodes)}")
    print(f"Total nodes: {gc.num_total}")
    
    print(f"\nConnectivity graph edges: {len(gc.connectivity_graph.edges())}")
    for i, j in gc.connectivity_graph.edges():
        print(f"  {i} -- {j}: {gc.all_nodes[i]} -- {gc.all_nodes[j]}")
    
    # Define a simple weight function (e.g., based on evolutionary likelihood)
    def weight_function(vec1, vec2, info1, info2):
        # Higher weight for transitions involving original nodes
        original_bonus = 0
        if info1['type'] == 'original':
            original_bonus += 1
        if info2['type'] == 'original':
            original_bonus += 1
        return 1.0 + original_bonus
    
    # Set weights and compute maximum spanning tree
    gc.set_edge_weights(weight_function)
    mst = gc.compute_maximum_weight_spanning_tree()
    
    print(f"\nMaximum Weight Spanning Tree:")
    tree_info = gc.get_tree_info()
    print(f"  Total weight: {tree_info['total_weight']}")
    print(f"  Tree edges:")
    for i, j, weight in tree_info['tree_edges']:
        print(f"    {i} -- {j} (weight: {weight})")
    
    # Test evolutionary path
    if len(original) >= 2:
        path = gc.get_evolutionary_path(0, len(original)-1)
        print(f"\nEvolutionary path from node 0 to node {len(original)-1}: {path}")


if __name__ == "__main__":
    test_graph_constructor()